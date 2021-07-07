import re

from common import *

#-----------------------------------------------------------------------------
# MDS IPSS Risk Classification
# Source: WCM Silver MPN Center 
#-----------------------------------------------------------------------------
'''
------------------------------------------------------------------------------
Risk Status   Single Abnormality  Double Abnormality  Complex Abnormality
------------------------------------------------------------------------------
Good          Normal                               
Good          -Y    
Good          del(5q)   
Good          del(20q)    
------------------------------------------------------------------------------
Intermediate  Other               Any 
------------------------------------------------------------------------------
Poor          -7    
Poor          del(7q)                                 count(abn)>2
------------------------------------------------------------------------------
'''

#-----------------------------------------------------------------------------

def good(k):
  res = [re_del_5q, re_del_20q]
  return is_normal(k) or (any([ searchall(res, k),missing_y(k) ])
                          and abn_count(k) == 1 )

#-----------------------------------------------------------------------------

def inter(k):
  '''Here, we need to check the value of good() but we cannot do the
  same with poor(),due to the "Complex Abnormality" rule; so we
  construct a list of some poor-specific tests that should be confirmed
  as false.
  (TODO confirm if this reasoning is correct.
  Technically, order of these calls matter.)'''
  poor_single_abnormality_res = [re_minus_7, re_del_7q]
  return (  (abn_count(k) == 1 
            and (not good(k))
            and (not searchall(poor_single_abnormality_res, k))) 
          or
            (abn_count(k) == 2))

#-----------------------------------------------------------------------------

def poor(k):
  res = [re_minus_7, re_del_7q]
  return any([
    searchall(res, k)
    ,abn_count(k) > 2 ])

#-----------------------------------------------------------------------------

def classify(k):
  k = k.strip().replace('\n', '')
  try:
    if poor(k): return 'poor'
    elif inter(k): return 'intermediate' 
    elif good(k): return 'good'
    return '?'
  except Exception, e:
    print str(e)
    return 'error'
 
#-----------------------------------------------------------------------------
# references and notes

'''
IPSS-R Cytogenetic Risk Groups:
  https://www.mds-foundation.org/ipss-r-calculator/

DIPSS Plus Score for Prognosis in Myelofibrosis:
  https://qxmd.com/calculate/calculator_315/dipss-plus-score-for-prognosis-in-myelofibrosis

'''

