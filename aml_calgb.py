import re
from common import *

#-----------------------------------------------------------------------------
# CALGB (Cancer and Leukemia Group B) Risk Classification
# Source: WCM Silver MPN Center 
#-----------------------------------------------------------------------------
'''
-------------------------------------------------------
Risk Status               CALGB
-------------------------------------------------------
Favorable                 t(8;21)
                          inv(16)/t(16;16)
-------------------------------------------------------
Intermediate              Normal
                          -Y
                          del(5q)
                          t(6;9)
                          t(6;11)
                          -7
                          del(7q)
                          +8
                          del(9q)
                          t(9;11)
                          +11
                          del(11q)
                          t(11;19)(q23;p.13.1)
                          +13
                          del(20q)
                          +21
-------------------------------------------------------
Unfavorable               inv(3)/t(3;3)
                          abn(12p)
                          count(abn)>2
-------------------------------------------------------
'''

#-----------------------------------------------------------------------------

def fav(k):
  '''
  Favorable               t(8;21)
                          inv(16)/t(16;16)
  '''
  res = [re_t_8_21, re_inv_16, re_t_16_16]
  return searchall(res, k)

#-----------------------------------------------------------------------------
def inter(k):
  '''
  -------------------------------------------------------
  Intermediate              Normal
                            -Y
                            del(5q)
                            t(6;9)
                            t(6;11)
                            -7
                            del(7q)
                            +8
                            del(9q)
                            t(9;11)
                            +11
                            del(11q)
                            t(11;19)(q23;p13.1)
                            +13
                            del(20q)
                            +21
  -------------------------------------------------------
  '''
  res = [re_del_5q, re_t_6_9, re_t_6_11, re_minus_7, re_del_7q, re_plus_8
        ,re_del_9q, re_t_9_11, re_plus_11, re_del_11q, re_plus_13
        ,re_del_20q, re_plus_21]
  return any([
     is_normal(k) 
    ,searchall(res, k)
    ,missing_y(k)
    # t(11;19)(q23;p.13.1)
    ,t_involving(k, 11, 'q', 23, 19, 'p', 13.1) ])

#-----------------------------------------------------------------------------
def unf(k):
  '''
  -------------------------------------------------------
  Unfavorable               inv(3)/t(3;3)
                            abn(12p)
                            count(abn)>2
  -------------------------------------------------------
  '''
  res = [re_inv_3, re_t_3_3]
  return any([
    searchall(res, k)
    ,abn_at(k, 12, 'p')
    ,abn_count(k) > 2 ])

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

