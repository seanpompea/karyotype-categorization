import sys
import csv
import time
import json
import codecs
import datetime
import multiprocessing
import collections

from kickshaws import *
from sqlsrvwrapper import *

import common as c

#-----------------------------------------------------------------------------

import aml_swog
import aml_calgb
import aml_mrc_1998
import aml_mrc_2010
import mds_ipss
import mds_5_group
import mf_dipss_plus

#-----------------------------------------------------------------------------

def slurp_json(filepath):
  #with codecs.open(filepath, 'r', 'utf_8') as f:
  with open(filepath, 'r') as f:
    return json.load(f)

def now():
  return str(datetime.datetime.now())

def make_ordered_map(keys, mp):
  '''Returns a map with key-value pairs ordered according to the keys
  argument; can filter unwanted keys by excluding the key name from keys.
  '''
  out = collections.OrderedDict()
  outkeys = [k for k in keys if mp.get(k,'*sentinel*') is not '*sentinel*']
  for k in outkeys:
    out[k] = mp[k]
  return out

# this includes extra columns used for testing/verification.
cols_expanded = [
 'mrn'
,'karyotype'
,'result_date'
,'aml_swog'
,'aml_calgb'
,'aml_mrc_1998'
,'aml_mrc_2010'
,'mds_ipss'
,'mds_5_group'
,'mf_dipss_plus'
,'is_normal'
,'abn_count'
,'full_numer_abn_tally'
,'trailing_mar_and_r_abn_tally'
]

#-----------------------------------------------------------------------------
# csv

def slurp_csv(fname):
  '''Reads a CSV file and returns a list of maps. Assumes CSV is UTF-8 and
  map's values will contain unicode-type strings.'''
  with open(fname, 'r') as f:
  #with codecs.open(fname, 'r', 'utf_8_sig') as f:
    reader = csv.DictReader(f)
    # Below, we convert each value from raw bytes into a unicode-type string
    # with assumption the bytes have UTF-8 encoding.
    return [{k:unicode(v, 'utf-8') for k, v in row.iteritems()}
            for row in reader]

def data2csv(fname, data):
  '''Takes a list of maps and writes to a CSV file. Attempts
  to encode map values that are unicode-type strings into UTF-8-encoded 
  bytes.'''
  data = map(lambda mp: make_ordered_map(cols_expanded, mp), data)
  print(data)
  with open(fname, 'w') as f:
    colnames = data[0].keys()
    wrtr = csv.DictWriter(f, fieldnames=colnames)
    wrtr.writeheader()
    # Below, we check if each value (in a key-value pair) is of type unicode,
    # and, if so, encode it into bytes using UTF-8 encoding.
    wrtr.writerows([{k:v.encode('utf-8') if type(v) is unicode else v
                    for k, v in row.iteritems()}
                    for row in data])
  return fname

#-----------------------------------------------------------------------------
# db

def retrieve_from_db(db_tag):
  '''Queries the MPN database and returns karyotype data.
  Uses the conneciton details defined in an external file: if 
  db_source_id is P02, will use enclave/p02-info.json.'''
  db_info = slurp_json('enclave/' + db_tag.lower() + '-info.json')
  qy = 'select mrn, result_date, karyotypes as karyotype '\
       'from dm_mpn..MPN_KARYOTYPES '\
       'where karyotypes is not null'
  return db_qy(db_info, qy)

#-----------------------------------------------------------------------------
# process

def apply_classification(data, classification_module):
  '''For each map in data, add key-value pair (or update if the kvp
  has been preloaded, as is the case when we use multiple procs)
  where the key is the classification_module's name and 
  value is result of calling classification_module.classify.
  Note: modifies data directly; does not create a new copy.
  Assumes each map has a key of 'karyotype'.
  Returns data for convenience.'''
  print 'Starting ' + classification_module.__name__
  for mp in data:
    strat_rslt = classification_module.classify(mp['karyotype'])
    mp[classification_module.__name__] = strat_rslt
  print 'Finished ' + classification_module.__name__ + ' ' + now()
  return data

def apply_f(data, f):
  for mp in data:
    rslt = f(mp['karyotype'].strip().replace('\n',''))
    mp[f.__name__] = rslt
  print 'Finished ' + f.__name__ + ' ' + now()
  return data

def process(data):
  '''Apply all risk classifications to all maps in data. 
  Modifies data directly; does not create a new copy.
  To skip a module, just comment it out below.
  Returns data for convenience.'''
  modules = [aml_swog
            ,aml_calgb
            ,aml_mrc_1998
            ,aml_mrc_2010
            ,mds_ipss
            ,mds_5_group
            ,mf_dipss_plus]
  map(lambda m: apply_classification(data, m), modules)
  funcs = [c.abn_count
          ,c.is_normal
          ,c.full_numer_abn_tally
          ,c.trailing_mar_and_r_abn_tally]
  map(lambda f: apply_f(data, f), funcs)
  return data

# TODO broken
def mp_process(data):
  '''Variation of the 'process' function; this uses multiple processes
  via the multiprocessing module.'''
  modules = [aml_swog
            ,aml_calgb
            ,aml_mrc_1998
            ,aml_mrc_2010
            ,mds_ipss
            ,mds_5_group
            ,mf_dipss_plus]
  # Preload key-value pairs so that the step of adding a new key
  # is prevented from happening when more than one proc is accessing
  # a map simultaneously.
  mgr = multiprocessing.Manager()
  managed_data = mgr.list([mgr.dict(mp) for mp in data])
  for module in modules:
    for mp in managed_data:
      mp[module.__name__] = ''
  # Set up procs and run.
  procs = []
  for module in modules:
    p = multiprocessing.Process(target=mp_apply_classification
                                ,args=(managed_data, module))
    procs.append(p)
  # Start them all
  for p in procs:
    p.start()
  # ...and wait for all to finish.
  for p in procs:
    p.join()
  #return data
  return managed_data.items()

#-----------------------------------------------------------------------------
# drivers

table1 = 'dm_mpn.dbo.karyotype_classifications'
table2 = 'dm_mpn.dbo.karyotype_classifications2'
table3 = 'dm_mpn.dbo.karyotype_classifications3'

def trunctable(db_tag, table_name):
  db_trunc_table(
    slurp_json('enclave/' + db_tag + '-info.json'), table_name)
  
def csv2csv(fname, use_mp=None):
  '''Read in a CSV which has, at minimum, a column named karyotype; process
  the classifications and output a new CSV with additional columns for
  the risk classifications.
  Returns the name of the new CSV.'''
  output_filename = fname + '-processed-' + str(int(time.time())) + '.csv'
  if use_mp:
    return data2csv(output_filename, mp_process(slurp_csv(fname)))
  else:
    return data2csv(output_filename, process(slurp_csv(fname)))

def csv2db(fname, db_tag, table_name, use_mp=None):
  '''Pass use_mp=True to use multiple procs with multiprocessing.'''
  data = slurp_csv(fname)
  print 'Loaded CSV ' + str(datetime.datetime.now())
  print 'Starting classifications... ' + str(datetime.datetime.now())
  if (use_mp):
    rslt = mp_process(data)
  else:
    rslt = process(data)
  print 'Processed classifications ' + str(datetime.datetime.now())
  print 'Starting to load into database... ' + str(datetime.datetime.now())
  tbl = 'dm_mpn..karyotype_classifications2'
  db_insert_many(slurp_json('enclave/' + db_tag + '-info.json'), table_name, rslt)
  print 'Finished. ' + str(datetime.datetime.now())

def db2db(src_db_tag, dest_db_tag, table_name, use_mp=None):
  print 'Database to query from: ' + src_db_tag
  print 'Database to put results into: ' + dest_db_tag
  data = retrieve_from_db(src_db_tag)
  print 'total rows: ' + str(len(data))
  print 'Starting classifications... ' + str(datetime.datetime.now())
  if (use_mp):
    rslt = mp_process(data)
  else:
    rslt = process(data)
  print 'Processed classifications ' + str(datetime.datetime.now())
  print 'Starting to load into database... ' + str(datetime.datetime.now())
  db_insert_many(slurp_json('enclave/' + dest_db_tag + '-info.json'), table_name, rslt)
  print 'Finished. ' + str(datetime.datetime.now())

def db2csv(db_tag, use_mp=None):
  data = retrieve_from_db(db_tag)
  print 'total rows: ' + str(len(data))
  print 'Starting classifications... ' + str(datetime.datetime.now())
  if (use_mp):
    rslt = mp_process(data)
  else:
    rslt = process(data)
  print 'Processed classifications ' + str(datetime.datetime.now())
  outname = now() + '.csv'
  return data2csv(outname, rslt)

#-----------------------------------------------------------------------------

def main():
  print sys.argv[1]
  csv2csv(sys.argv[1])

if __name__ == '__main__': main()

#-----------------------------------------------------------------------------
'''

A typical load workflow might look like this:

>>> import bulk as b
>>> b.trunctable('p06', 'dm_mpn..karyotype_classifications')
>>> b.db2db('p06', 'p04', 'dm_mpn..karyotype_classifications')

'''


