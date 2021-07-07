import functools, itertools, operator
import re
from common import *

#------------------------------------------------------------------------------
# AML SWOG (Southwest Oncology Group) Risk Classification
# Source: WCM Silver MPN Center
#------------------------------------------------------------------------------
'''
-------------------------------------------------
Risk Status       Abnormalities
-------------------------------------------------
Favorable         t(15;17) 
                  t(8;21) 
                  inv(16)/t(16;16)/del(16q) 
-------------------------------------------------
Intermediate      Normal
                  +8
                  +6
                  -Y
                  del(12p)
-------------------------------------------------
Unfavorable       abn(3q)
                  del(5q)/-5
                  -7/del(7q)
                  t(6;9)
                  t(9;22)
                  9q
                  11q
                  20q
                  21q
                  17p
                  count(abn)>2
-------------------------------------------------
Unknown           Other
-------------------------------------------------
'''

#------------------------------------------------------------------------------

def fav(k):
  res = [re_t_15_17, re_t_8_21, re_inv_16, re_t_16_16, re_del_16q]
  return searchall(res, k)

#------------------------------------------------------------------------------

def inter(k):
  res = [re_plus_8, re_plus_6, re_del_12p]
  return any([
    is_normal(k)
    ,missing_y(k)
    ,searchall(res, k) ])

#------------------------------------------------------------------------------

def unf(k):
  '''
  Unfavorable       abn(3q)
                    del(5q)/-5
                    -7/del(7q)
                    t(6;9)
                    t(9;22)
                    9q
                    11q
                    20q
                    21q
                    17p
                    count(abn)>2
  '''
  res = [re_minus_5, re_del_5q, re_minus_7, re_del_7q, re_t_6_9, re_t_9_22]
  return any([
    searchall(res, k)
    ,abn_at(k, 3, 'q')
    ,abn_at(k, 9, 'q')
    ,abn_at(k, 11, 'q')
    ,abn_at(k, 20, 'q')
    ,abn_at(k, 21, 'q')
    ,abn_at(k, 17, 'p')
    ,(abn_count(k) > 2) ])

#------------------------------------------------------------------------------

def classify(k):
  k = k.strip().replace('\n', '')
  try:
    if unf(k): return 'unfavorable'
    elif inter(k): return 'intermediate'
    elif fav(k): return 'favorable'  
    else: return 'unknown'
  except Exception, e:
    print str(e)
    return 'error'

