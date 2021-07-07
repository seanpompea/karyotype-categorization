import re

from common import *

#-----------------------------------------------------------------------------
# MF DIPSS+ Risk Classification
# Source: WCM Silver MPN Center 
#-----------------------------------------------------------------------------
'''
------------------------------------------------------------------------------
Risk Status             DIPSS+
------------------------------------------------------------------------------
Unfavorable             +8
                        -7/7q-
                        i(17q)
                        -5/5q-
                        12p-
                        inv(3)
                        11q23
------------------------------------------------------------------------------
Favorable               Other
------------------------------------------------------------------------------
'''

#-----------------------------------------------------------------------------

def unf(k):
  res = [re_plus_8, re_minus_7, re_del_7q, re_i_17q, re_minus_5
        ,re_del_5q, re_del_12p, re_inv_3]
  return any([searchall(res, k), abn_involving_11q23(k)])

#-----------------------------------------------------------------------------

def classify(k):
  k = k.strip().replace('\n', '')
  try:
    if unf(k): return 'unfavorable'
    else: return 'favorable'
  except Exception, e:
    print str(e)
    return 'error'

