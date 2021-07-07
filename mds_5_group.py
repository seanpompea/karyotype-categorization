import re

from common import *

#-----------------------------------------------------------------------------
# MDS 5-Group Risk Classification
# Source: WCM Silver MPN Center 
#-----------------------------------------------------------------------------
'''
------------------------------------------------------------------------------
Risk Status   Single Abnormality    Double Abnormality  Complex Abnormality
------------------------------------------------------------------------------
Very good     -Y    
Very good     del(11q)    
------------------------------------------------------------------------------
Good          Normal                del(5q) 
Good          del(5q)   
Good          del(20q)    
Good          del(12p)    
------------------------------------------------------------------------------
Intermediate  del(7q)               Any other 
Intermediate  +8    
Intermediate  i(17q)    
Intermediate  +19   
Intermediate  Any other   
------------------------------------------------------------------------------
Poor          -7    
Poor          inv(3)/t(3q)/del(3q)  -7/del(7q)          count(abn)=3
------------------------------------------------------------------------------
Very poor                                               count(abn)>3
------------------------------------------------------------------------------
'''

#-----------------------------------------------------------------------------

def vg(k):
  return (any([
      missing_y(k)
      ,re_del_11q.search(k)])
    and abn_count(k) == 1)

#-----------------------------------------------------------------------------

def good(k):
  res = [re_del_5q, re_del_20q, re_del_12p]
  return (
    is_normal(k)
    or
    (searchall(res, k) and abn_count(k) == 1)
    or
    (re_del_5q.search(k) and abn_count(k) == 2) )

#-----------------------------------------------------------------------------

def inter(k):
  res = [re_del_7q, re_plus_8, re_plus_19, re_i_17q]
  return ( 
    (any([searchall(res, k), abn_at(k, 17, 'q')]) and abn_count(k) == 1)
    or
    (abn_count(k) == 2) )

#-----------------------------------------------------------------------------

def poor(k):
  res_single_abn = [re_minus_7, re_inv_3, re_del_3q] 
  res_dbl_abn = [re_minus_7, re_del_7q]
  return ( 
    (any([searchall(res_single_abn, k), t_involving_3q(k)]) and abn_count(k) == 1)
    or
    (searchall(res_dbl_abn, k) and abn_count(k) ==2)
    or
    (abn_count(k) == 3) )

#-----------------------------------------------------------------------------

def vp(k):
  return abn_count(k) > 3

#-----------------------------------------------------------------------------

def classify(k):
  k = k.strip().replace('\n', '')
  try:
    if vp(k): return 'very poor'
    elif poor(k): return 'poor'
    elif inter(k): return 'intermediate'
    elif good(k): return 'good'
    elif vg(k): return 'very good'
    return 'intermediate'
  except Exception, e:
    print str(e)
    return 'error'

