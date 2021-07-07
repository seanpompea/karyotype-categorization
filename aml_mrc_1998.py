import re
from common import *

#-----------------------------------------------------------------------------
# AML MRC 1998 Risk Classification
# Source: WCM Silver MPN Center 
#-----------------------------------------------------------------------------
'''
-----------------------------------------------------------------
Risk Status                 MRC 1998      
-----------------------------------------------------------------
Favorable                   t(15;17)
                            t(8;21)
                            inv(16)/t(16;16)/del(16q)
-----------------------------------------------------------------
Intermediate                Normal
                            11q23 abn
                            +8
                            del(9q)
                            del(7q)
                            +21
                            +22
-----------------------------------------------------------------
Unfavorable                 abn(3q)
                            del(5q)/-5
                            -7
                            count(abn)>4
-----------------------------------------------------------------
Note: table in article belo says 'all others' should be classified
as intermediate; this is not specified by Silver however.
Orozco et al: http://www.cancernetwork.com/acute-myeloid-leukemia/unfavorable-complex-and-monosomal-karyotypes-most-challenging-forms-acute-myeloid-leukemia
'''


#-----------------------------------------------------------------------------

def fav(k):
  res = [re_t_15_17, re_t_8_21, re_inv_16, re_t_16_16, re_del_16q]
  return searchall(res, k)

#-----------------------------------------------------------------------------

def inter(k):
  '''
  ----------------------------------------------
  Intermediate                Normal
                              11q23 abn
                              +8
                              del(9q)
                              del(7q)
                              +21
                              +22
  ----------------------------------------------
  '''
  simple_res = [re_plus_8, re_del_9q, re_del_7q, re_plus_21
               ,re_plus_22]
  return any([ 
    is_normal(k)
    ,abn_involving_11q23(k)
    ,searchall(simple_res, k) ])

#-----------------------------------------------------------------------------

def unf(k):
  '''
  ----------------------------------------------
  Unfavorable                 abn(3q)
                              del(5q)/-5
                              -7
                              count(abn)>4
  ----------------------------------------------
  '''
  simple_res = [re_minus_5, re_del_5q, re_minus_7]
  return any([ 
    abn_at(k, 3, 'q')
    ,searchall(simple_res, k)
    ,abn_count(k) > 4 ])

#-----------------------------------------------------------------------------

def classify(k):
  k = k.strip().replace('\n', '')
  try:
    if unf(k): return 'unfavorable'
    elif inter(k): return 'intermediate'
    elif fav(k): return 'favorable'  
    return '?'
  except Exception, e:
    print str(e)
    return 'error'
 


