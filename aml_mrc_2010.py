import re

from common import *

#-----------------------------------------------------------------------------
# AML MRC 2010 Risk Classification
# Source: WCM Silver MPN Center 
#-----------------------------------------------------------------------------
'''
------------------------------------------------------------------------------
Risk Status             MRC 2010      
------------------------------------------------------------------------------
Favorable               t(15;17)(q22;q21)
                        t(8;21)(q22;q22)
                        inv(16)(p13q22)/t(16;16)(p13;q22)
------------------------------------------------------------------------------
Intermediate                
  If does not fall into favorable or unfavorable, classify as intermediate.
  See logic in classify function.
------------------------------------------------------------------------------
Unfavorable             abn(3q) [excluding t(3;5)(q21~25;q31~35)]
                        inv(3)(q21q26)/t(3;3)(q21;q26)
                        add(5q)
                        del(5q)
                        -5
                        add(7q)/del(7q)
                        -7
                        t(6;11)(q27;q23)
                        t(10;11)(p11~13;q23)
                        t(11q23) [excluding t(9;11)(p21~22;q23) and t(11;19)(q23;p13)]
                        t(9;22)(q34;q11)
                        -17/abn(17p)
                        count(abn)>3
------------------------------------------------------------------------------
'''

def fav(k):
  return any([
    re_inv_16_p13q22.search(k)
    # t(15;17)(q22;q21)
    ,t_involving(k, 15, 'q', 22, 17, 'q', 21)
    # t(8;21)(q22;q21)
    ,t_involving(k, 8, 'q', 22, 21, 'q', 22)
    # t(16;16)(p13;q22)
    ,t_involving(k, 16, 'p', 13, 16, 'q', 22) ])

#-----------------------------------------------------------------------------

# Note:
# this module doesn't have an 'intermediate' function.

#-----------------------------------------------------------------------------
# unfavorable

def aml2010_t_exclusion_1(k):
  '''This returns true if the kar matches t(3;5)(q21~25;q31~35)]'''
  tests = []
  for i in range(21,25):
    for j in range(31,36):
      tests.append(t_involving(k, 3, 'q', i, 5, 'q', j))
  return any(tests)

def aml2010_t_exclusion_2(k):
  '''This returns true if the kar matches
  t(9;11)(p21~22;q23) or t(11;19)(q23;p13)]'''
  return any([
    t_involving(k, 9, 'p', 21, 11, 'q', 23)
    ,t_involving(k, 9, 'p', 22, 11, 'q', 23)
    ,t_involving(k, 11, 'q', 23, 19, 'p', 13)])

def unf(k):
  '''
  ----------------------------------------------
  Unfavorable       abn(3q) [excluding t(3;5)(q21~25;q31~35)]
                    inv(3)(q21q26)/t(3;3)(q21;q26)
                    add(5q)
                    del(5q)
                    -5
                    add(7q)/del(7q)
                    -7
                    t(6;11)(q27;q23)
                    t(10;11)(p11~13;q23)
                    t(11q23) [excluding t(9;11)(p21~22;q23) and t(11;19)(q23;p13)]
                    t(9;22)(q34;q11)
                    -17/abn(17p)
                    count(abn)>3
  ----------------------------------------------
  '''
  return any([
    # abn(3q) [excluding t(3;5)(q21~25;q31~35)]
    (abn_at(k, 3, 'q') and not aml2010_t_exclusion_1(k))

    # inv(3)(q21q26)
    ,re_inv_3_q21q26.search(k)

    # t(3;3)(q21;q26)
    ,t_involving(k, 3, 'q', 21, 3, 'q', 26) 

    # add(5q)
    ,re_add_5q.search(k) 

    # del(5q)
    ,re_del_5q.search(k)

    # -5
    ,re_minus_5.search(k)

    # add(7q)
    ,re_add_7q.search(k)

    # del(7q)
    ,re_del_7q.search(k)

    # -7
    ,re_minus_7.search(k)

    # t(6;11)(q27;q23)
    ,t_involving(k, 6, 'q', 27, 11, 'q', 23)

    # t(10;11)(p11~13;q23)
    ,t_involving(k, 10, 'p', 11, 11, 'q', 23)
    ,t_involving(k, 10, 'p', 12, 11, 'q', 23)
    ,t_involving(k, 10, 'p', 13, 11, 'q', 23)

    # t(11q23) [excluding t(9;11)(p21~22;q23) and t(11;19)(q23;p13)]
    ,(t_involving_11q23(k) and not aml2010_t_exclusion_2(k))

    ,(abn_count(k) > 3) ])

#-----------------------------------------------------------------------------

def classify(k):
  k = k.strip().replace('\n', '')
  try:
    if unf(k): return 'unfavorable'
    elif fav(k): return 'favorable'  
    return 'intermediate'
  except Exception, e:
    print str(e)
    return 'error'
 



