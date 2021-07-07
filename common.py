from functools import *
from itertools import *
from operator import *

import re

#-----------------------------------------------------------------------------
# conventions in the code
'''

're' vs 'plain'
--------------
In variable/function names, "re_" is used to indicate a string/slug 
or regex use versus "plain_" to indicate non-regex use. 

'setup__' functions
-------------------
A pattern is used here where a function establishes some compiled
regex objects (perhaps in a cache) locally and then creates a closure by 
returning an inner function.
The "outer"/callable function is then instantiated immediately after.
This is for performance and also to keep regexes close syntactically
to the location where they are used.

'''

#----------------------------------------------------------------------------- 

def searchall(res, k):
  '''Takes a list of compiled gegex objects, calls search on them all
  and returns true if any find a match.'''
  return any(map(lambda x: x.search(k), res)) 

def findall(res, k):
  '''Similar to above but uses 'findall' instaed of 'search'.'''
  return any(map(lambda x: x.findall(k), res)) 

#-----------------------------------------------------------------------------
# useful regex slugs

re_slug_before = r'(?:[\s\,\)\+])'

re_slug_after =  r'(?:[\(\s\,\[\.]|$)'
  # This is more liberal than 'slug_after_all' below.

#TODO might have other stuff ('e.g., '.ish') trailing the count, so we
# need to accomodate that.
re_slug_after_all = r'(?:[\s\,\[\.]|$)'
  # i.e., nearing the end of the karyotype string; no more events, only 
  # cell count and/or additional trailing modifiers etc. 

#-----------------------------------------------------------------------------
# numerical events (see AGT Section 8.4 "Numerical Events" p378)

def setup__full_numer_abn_tally(): # arg will be just the karyotype string
  '''Tally of 'full' numerical abnormalities (E.g., '+21'). This function
  does not attempt to find other types of numerical abnormalities involving
  a non-complete/partial homologue etc.

  Note: for accuracy, only use 'abn_count' function in production cases..

  This fn does not find 'partial' scenarios -- eg partial trisomy, etc.

  Finds explicit full/entire extra/missing chromosomes (sex or autosome), 
  along with implicit extra/missing sex chromosomes (implicit example: 45,X).
  '''
  re1 = re.compile(r'\-(?:X|Y)')
  re2 = re.compile(r'\d{2}\,(?:X|Y)\b')
  re3 = re.compile(r'\d{2}\,(?:X|Y){3,}')
  re4 = re.compile(r'[-+]\d+' + re_slug_after)
  # The above will be cached and not rebuilt; we create a closure by
  # returning an inner functio.
  def full_numer_abn_tally(k): # Use full name so __name__ is ok later.
    '''Discover *explicit* missing sex chromosome.
      >>> re.findall('\-(?:X|Y)', '46,XX') []
      >>> re.findall('\-(?:X|Y)', '46,X,-Y')
      ['-Y']
    '''
    c1 = re.findall(re1, k)
    '''Discover *implicit* missing sex chromosome.
      >>> re.findall(r'\d{2}\,(?:X|Y)\b', '45,Y')
      ['45,Y']
      >>> re.findall(r'\d{2}\,(?:X|Y)\b', '45,XX')
      []
      >>> re.findall(r'\d{2}\,(?:X|Y)\b', '45,XY')
      []
    '''
    c2 = []
    if not c1: # Only look if not explicit.
      c2 = re.findall(re2, k)
    '''Discover additional sex chromosome(s).
      >>> re.findall(r'\d{2}\,(?:X|Y){3,}', '45,XXXY')
      ['45,XXXY']
    '''
    c3 = re.findall(re3, k)
    '''Discover missing or extra autosome homologues.
      >>> re.findall('[-+]\d+', '48,XX,+18,+21')
      ['+18', '+21']
    '''
    c4 = re.findall(re4, k)
    return sum(map(len, [c1,c2,c3,c4]))
  return full_numer_abn_tally
# instantiate function for use.
full_numer_abn_tally = setup__full_numer_abn_tally()

#------------------------------------------------------------------------------

def setup__trailing_mar_and_r_abn_tally(): # arg is just kar string.
  '''Look for 'mar' and 'r' at the end of a karyotype. Both symbols might be 
  subscripted if there are clonally distinct instances, and can be preceded
  by a number (indicated count) and/or a '+' sign (indicating supernumerary). 
  This fn does *not* look for rings of known origin, which has the 
  typical r(N)... format.
  '''
  re1 = re.compile(re_slug_before + r'\+?\d?mar')
  re2 = re.compile(r'[^\w]\+?\d?r')
  def trailing_mar_and_r_abn_tally(k):
    '''Step one: find 'mar' instances.
      >>> k1
      '47,XY,+3,-7,-20,+der(?)t(?;1)(?;p11.2),+mar,inc[20]'
      >>> k2
      '47,XY,+3,-7,-20,+der(?)t(?;1)(?;p11.2),mar,inc[20] '
      >>> k3
      '47,XY,+3,-7,-20,+der(?)t(?;1)(?;p11.2),mar1,mar2,inc[20] '
      >>> re.findall(r'\+?\d?mar', k1)
      ['+mar']
      >>> re.findall(r'\+?\d?mar', k2)
      ['mar']
      >>> re.findall(r'\+?\d?mar', k3)
      ['mar', 'mar']
    '''
    #TODO do we need the before slug?
    c1 = re.findall(re1, k)
    '''Step two: find 'r' instances near end of kar -- these are rings
    of unknown origin, (as opposed to r(N) which provides known
    chromosome number)
    Note that this one also as a  'der' symbol which we don't want to
    inadvertently count here.
      >>> k5
      '47,XY,+3,-7,-20,+der(?)t(?;1)(?;p11.2),+2r1,+r2,inc[20] '
      >>> re.findall(r'[^\w]\+?\d?r', k5)
      [',+2r', ',+r']
    '''
    c2 = re.findall(re2, k)
    return sum(map(len, [c1, c2]))
  return trailing_mar_and_r_abn_tally
trailing_mar_and_r_abn_tally = setup__trailing_mar_and_r_abn_tally()

#------------------------------------------------------------------------------

def setup__inc_abn_tally(): # arg is just kar string.
  '''Abnormality tally for 'inc' symbol.
  Returns return 1 or 0 ( unless composite). 'inc' or '+inc' indicates 
  additional abnormality material that could not be identified.'''
  re1 = re.compile(r'\+?inc' + re_slug_after_all)
  def f(k):
    return len(re.findall(re1, k))
  return f
inc_abn_tally = setup__inc_abn_tally()

#------------------------------------------------------------------------------

def setup__dmin_abn_tally():  # arg is just kar string.
  '''Abnormality tally for 'dmin' symbol.
  Returns return 1 or 0 (unless composite). dmin == double minute.
  Examples:
    46,XY,5dmin[20]
    46,XY,5~ 22dmin[20]  <-- Not sure if AGT intends the space; leaving as-is
                             though.
    -- two above From AGT p402.
    46,XX,2~15dmin[20] (intra; AGT 424)
  '''
  r = re.compile(r'\~?\s?\d{1,2}dmin' + re_slug_after_all)
  def f(k):
    return len(re.findall(r'\~?\s?\d{1,2}dmin' + re_slug_after_all, k))
  return f
dmin_abn_tally = setup__dmin_abn_tally()

#------------------------------------------------------------------------------
# structural abnormalities
'''
Event symbols are broken down here into several categories:
  o intrachromosomal
  o interchromosomal
  o janus (can be intra or inter)
  o derivative
'''

intra_events = {
  # Primary source: AGT Section 8.5 "Structural events" p380  
  # Going with AGT, ,the below contains "most structural events. Because of
  # their potential complexity, derivatives and symbols of uncertainty will be
  # discussed in their own sections." (AGT p380)
   'add':   'addition'                # AGT 8.7.2 p 399  *later section
  ,'del':   'deletion'                # AGT 8.5.1 p380
  ,'dup':   'duplication'             # AGT 8.5.3 p383
  ,'trp':   'triplication'            # AGT 8.5.3 p384
  ,'qdp':   'quadruplication'         # AGT 8.5.3 p 384
  ,'fis':   'cenral fission'          # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2693553/
  ,'fra':   'rare fragile site'       # AGT 8.8.4 p 404; ex: 46,Y,fra(X)(q27.3)
  ,'hsr':   'Homogeneously staining region' # AGT 8.7.2 p 399 *later section 
  ,'inv':   'inversion'               # AGt 8.5.5 p 385
                                      # Unusual example w/ sex chrom: inv(Y)(p11.2q11.2)  (AGT p404)
  ,'inv dup': 'inverted duplication'  #  CA p400 eg 46,XY, inv dup(1)(q32q21)[20] 
                                      # TODO: dcsn: might not need to do this separately
  ,'i':     'isochromosome'           # AGT 8.5.6 p 386
  ,'idic':  'isodicentric'            # AGT 8.5.7 p 387 
  ,'psu idic':'pseudoisodicentric'    # AGT 387
  ,'psudic': 'pseudoisodicentric'     # ... alternate form per Cooper 160
                                      # AGT 8.5.8 consitutional annotators; see below.
  #,'rec':   'recombinant'             # AGT 8.5.9 p 389
  # Analysis of WCM data: no usage of 'rec'
  ,'r':     'ring of known centric origin'
                                      # AGT 8.5.10 p 390
                                      # See also: 'dic r' and 'trc c' below.
  # AGT 8.7.6 p 400  -- ring of unknown origin: +r occur by itself
  #   eg: 47,XY,-18,+21,+r
  # See below, +r below in 'uncertainty'
  ,'tas':   'telomeric association'   #  AGT 8.5.11 p 391; multi
  ,'upd':   'uniparental disomy'      # AGT 8.5.13 p 394; eg 46,XY,upd(15)pat
}

inter_events = {
   'dic':   'dicentric'               # AGT 8.5.2 p381
  ,'psu dic': 'pseudodicentric'       # AGT 8.5.2 p383
  ,'dic r' : 'dicentic ring'          # AGT -- multi: two chromosomes involved
  ,'trc r' : 'tricentic ring'         # AGT -- multi: 3 chromosomes involved
  ,'t':     'translocation'           # AGT 8.5.12 p 391; multi
}

janus_events = {
  # Events that might be intra or inter. Eesiest to 
  # group these separately.
  'ins':   'insertion'               # AGT 8.5.4 p 385 
}

derivative_events = {
   'der': 'derivative'
  ,'ider': 'isoderivative'
  #,'neo': 'neoderivative' # turn off for now. not common.
}

def der_preceding_qualifying_event(k):
  '''Returns true if k contains a der() which is then followed
  by a qualifying event. (Other times, der() stands alone, and this
  function should return False.)
  E.g.  2-event derivative:
          46,XY,der(2)t(2;5)(p23;q35)inv(2)(q21q31)
  (AGT p425)'''
  #TODO
  pass

re_event_prefixes = ['\)', '\+', ',', ' ']
  # + a plus indicates a kind of (partial) trisomy
  # We escape special chars having special meaning in a regex context.

def gen_abn_struct_event_re_slugs():
  slugs = []
  x = re_event_prefixes
  y = {}
  y.update(intra_events)
  y.update(inter_events)
  y.update(janus_events)
  y.update(derivative_events)
  for a in x:
    for b in y:
      slugs.append(a + b + r'\(')
  return slugs

# Use these slugs when constructing a regex.
struct_event_re_slugs = gen_abn_struct_event_re_slugs()

plain_event_prefixes = [')', '+', ',', ' ']
  # not for regex use; rather just plain slug/substring use.
  # + a plus indicates a kind of (partial) trisomy
  
def gen_abn_struct_event_plain_slugs():
  slugs = []
  x = plain_event_prefixes
  y = {}
  y.update(intra_events)
  y.update(inter_events)
  y.update(janus_events)
  y.update(derivative_events)
  for a in x:
    for b in y:
      slugs.append(a + b + r'(') # no escaping the (
  return slugs

# Use these when you just want to search a sring
# (e.g., using the find() function)
struct_event_plain_slugs = gen_abn_struct_event_plain_slugs()

def abn_count(k):
  # First look for structural event symbols.
  slugs = struct_event_plain_slugs
  c = 0
  for slug in slugs:
    #if k.find(slug) != -1:
    #  c += 1
    c+= k.count(slug)
  # Now add this up with the results of the other sub-tally funcs (which are
  # above.
  # TODO subtract 1 where der where has associated. See func above.
  return (c 
          + full_numer_abn_tally(k) 
          + trailing_mar_and_r_abn_tally(k)
          + inc_abn_tally(k)
          + dmin_abn_tally(k)
        )

#-----------------------------------------------------------------------------
# positional / structural / band involved

def setup__abn_at(): # args: kar, chrom, arm
  # TODO is a trisomy/monos of the chromosome relevant?
  '''
  Returns true if the karyotype has an abnormality involving the
  chromosome and arm.

  Notes:

  o "When more than one independent, structural event has occurred, the
    karyotype event must be written as a derivative" (AGT 8.6.4 p 395).
  '''
  #=================================================
  cache = {}
  def get_re(slug, chrom, arm, txt):
    '''This function facilitates compilation, caching, and retrieval of 
    regex objects dynamically just for the abn_at function.'''
    ky = slug + chrom + arm + txt
    if cache.get(ky):
      return cache[ky]
    else:
      r = None
      if txt == 'intra1':
        r = re.compile(slug + chrom + r'\)\(' + arm)
      elif txt == 'intra2':
        r = re.compile(slug + chrom + r'\)\([p|q]\d\d' + arm)
      elif txt == 'inter1':
        r = re.compile(slug + chrom + r';\s?\d{1,2}' + r'\)\(' + arm)
      elif txt == 'inter2':
        r = re.compile(slug + r'\d{1,2}\;\s?' + chrom + r'\)\([p|q]\d\d;\s?' + arm)
      else:
        raise Exception('unknown txt passed to get_re: ' + txt) 
      cache[ky] = r
      return r
  #=================================================
  def f(k, chrom, arm): # this func to be returned 
    chrom = str(chrom)
    rslts = []
    # TODO review this commented-out code.
    '''
    Remove the below logic for now; assumption had been that a numerical
    event is counted here, but it's not enough to identify this, 
    eg:
      '47,XY,+7,del(7)(p14)pat*'
    here, there is an extra homologue with a structural abnormality; the
    q arm is not invovled so if we only check for "+7' we might mistakenly
    return true. 

    c1 = []
    c2 = []
    c3 = []
    if chrom == 'Y' or chrom == 'X':
      # 1) Is it a sex chromosome and is it explicitly missing?
      c1 = re.findall('\-(?:' + chrom + ')', k)
      # 2) Sex chromosome and implicitly missing?
      c2 = re.findall(r'\d{2}\,(?:' + chrom + r')\b', k)
      # Additional sex chrom?
      if chrom == 'X':
        # Look for 3 or more
        c3 = re.findall(r'\d{2}\,(?:' + chrom + '){3,}', k)
      if chrom == 'Y':
        # Look for 2 or more
        c3 = re.findall(r'\d{2}\,(?:' + chrom + '){2,}', k)

    # 3) Is it a monosomy/trisomy? (Chromosom is entirely missing?)
    c4 = re.findall('[-+]' + chrom + slug_after, k)
    '''
    # 4) Look for the chrom/arm combination in tandem with a structural
    #    event symbol.
    tests = []
    # TODO, for relevant tests below,  
    # account for when breakpoints specify subband, eg q11.2.
    for slug in struct_event_re_slugs:
      # Recall that prepped slug will be ')inv(' or ',inv(' etc.
      # intra, arm of interest is immediately inside paren.
      tests.append(re.search(get_re(slug, chrom, arm, 'intra1'), k))
      # intra, arm of interest is part of a 2nd breakpoint
      tests.append(re.search(get_re(slug, chrom, arm, 'intra2'), k))
      # inter, breakpoint of interest is on first chrom
      tests.append(re.search(get_re(slug, chrom, arm, 'inter1'), k))
      # inter, arm of interest is part of 2nd chromosome
      tests.append(re.search(get_re(slug, chrom, arm, 'inter2'), k))
    #return (sum(map(len, [c1,c2,c3,c4])) or any(tests))
    return any(tests)
  return f
abn_at = setup__abn_at()

#------------------------------------------------------------------------------

def gen_t_re_slugs():
  slugs = []
  for a in re_event_prefixes:
    slugs.append(a + 't' + '\(')
  return slugs

t_re_slugs = gen_t_re_slugs()

def gen_inv_re_slugs():
  slugs = []
  for a in re_event_prefixes:
    slugs.append(a + 'inv' + '\(')
  return slugs

inv_re_slugs = gen_inv_re_slugs()

# TODO refactoring: next funcs are similar

def t_involving(k, ch1, arm1, bp1, ch2, arm2, bp2):
  '''Returns true if has translocation at chrome/arm/band combination.
    TODO  multichromosomal translocations eg
     46,XY,t(3;12;8;21)(q31;p13;q22;q22)  (AGT 393)
     46,XY,t(5;12;3;21)(q31;p13;p11;q22)  (AGT 393)
  '''
  ch1 = str(ch1)
  bp1 = str(bp1)
  ch2 = str(ch2)
  bp2 = str(bp2)
  tests = []
  for slug in t_re_slugs:
    r1 = (slug + ch1 + ';\s?' + ch2 + r'\)\s?\(' 
            + arm1 + bp1 + '\.?\d{0,2}?' + ';\s?' 
            + arm2 + bp2 + '\.?\d{0,2}\)')  
    tests.append(re.search(r1, k))
  return any(tests)

def setup__abn_involving_11q23(): # arg is kar string.
  #=================================================
  cache = {}
  def get_re(slug, txt):
    '''This function facilitates compilation, caching, and retrieval of 
    regex objects dynamically just for abn_involving_11q23.'''
    ky = slug + txt
    if cache.get(ky):
      return cache[ky]
    else:
      r = None
      if txt == 'intra1':
        r = re.compile(slug + r'11\)\(q23')
      elif txt == 'intra2':
        r = re.compile(slug + '11' + r'\)\([p|q]\d\d\.?\d{0,2}?' + 'q23')
      elif txt == 'inter1':
        r = re.compile(slug + '11;\s?' + r'\d{1,2}\)\(q23')
      elif txt == 'inter2':
        r = re.compile(slug + r'\d{1,2}\;\s?' + '11' 
                       + r'\)\([p|q]\d\d\.?\d{0,2}?;\s?' + 'q23')
      else:
        raise Exception('unknown txt passed to get_re: ' + txt)
      cache[ky] = r
      return r
  #=================================================
  def f(k):
    tests = []
    for slug in struct_event_re_slugs:
      # Recall that prepped slug will be ')inv(' or ',inv(' etc.
      # intra, arm of interest is immediately inside paren.
      tests.append(re.search(get_re(slug, 'intra1'), k))
      # intra, arm of interest is part of a 2nd breakpoint
      tests.append(re.search(get_re(slug, 'intra2'), k))
      # inter, breakpoint of interest is on first chrom
      tests.append(re.search(get_re(slug, 'inter1'), k))
      # inter, arm of interest is part of 2nd chromosome
      tests.append(re.search(get_re(slug, 'inter2'), k))
    return any(tests)
  return f
abn_involving_11q23 = setup__abn_involving_11q23()

def t_involving_11q23(k):
  '''This is similar to 'abn_involving_11q23 but just for translocation.'''
  tests = []
  for slug in t_re_slugs:
    # inter, breakpoint of interest is on first chrom
    tests.append(re.search(slug + '11;\s?' + r'\d{1,2}\)\(q23', k))
    # inter, arm of interest is part of 2nd chromosome
    tests.append(re.search(slug + r'\d{1,2}\;\s?' + '11' 
                           + r'\)\([p|q]\d\d\.?\d{0,2}?;\s?' + 'q23', k))
  return any(tests)

def t_involving_3q(k):
  tests = []
  for slug in t_re_slugs:
    # inter, breakpoint of interest is on first chrom
    tests.append(re.search(slug + '3;\s?' + r'\d{1,2}\)\(q', k))
    # inter, arm of interest is part of 2nd chromosome
    tests.append(re.search(slug + r'\d{1,2}\;\s?' + '3' 
                           + r'\)\([p|q]\d\d\.?\d{0,2}?;\s?' + 'q', k))
  return any(tests)

#-----------------------------------------------------------------------------
# convenience tests

res_normal = [re.compile(r'^46,XX(?:\[\d+\])?$')
             ,re.compile(r'^46,XY(?:\[\d+\])?$')]
  #TODO when we look into fish, review this alternative: 
  #   r'^46,XY(?:\[\d+\])(?:$| nuc ish|\.ish )'

def is_normal(k):
  return searchall(res_normal, k)

res_missing_y = [
  # -Y explicit --  e.g., '45,X,-Y'
  re.compile(r',-Y') 
  # -Y implicit -- e.g., '45,X'
  ,re.compile(r',X' + re_slug_after) ]

def missing_y(k):
  return searchall(res_missing_y, k)

#-----------------------------------------------------------------------------
# regex objects 
# prepared ahead of time for performance (among other reasons)

# +6
re_plus_6 = re.compile(r',\s?\+6' + re_slug_after)

# +8
re_plus_8 = re.compile(r',\s?\+8' + re_slug_after)

# +11
re_plus_11 = re.compile(r',\s?\+11' + re_slug_after)

# +13
re_plus_13 = re.compile(r',\s?\+13' + re_slug_after)

# +19
re_plus_19 = re.compile(r',\s?\+19' + re_slug_after)

# +21
re_plus_21 = re.compile(r',\s?\+21' + re_slug_after)

# +22
re_plus_22 = re.compile(r',\s?\+22' + re_slug_after)

# -5
re_minus_5 = re.compile(r',\s?-5' + re_slug_after)

# -7
re_minus_7 = re.compile(r',\s?-7' + re_slug_after)

# add(5q)
re_add_5q = re.compile(r'add\(5\)\(q')

# add(7q)
re_add_7q = re.compile(r'add\(7\)\(q')

# del(3q)
re_del_3q = re.compile(r'del\(3\)\(q')

# del(5q)
re_del_5q = re.compile(r'del\(5\)\(q')

# del(7q)
re_del_7q = re.compile(r'del\(7\)\(q')

# del(9q)
re_del_9q = re.compile(r'del\(9\)\(q')

# del(11q)
re_del_11q = re.compile(r'del\(11\)\(q')

# del(12p)
re_del_12p = re.compile(r'del\(12\)\(p')

# del(16q)
re_del_16q = re.compile(r'del\(16\)\(q')

# del(20q)
re_del_20q = re.compile(r'del\(20\)\(q')

# del(20q)
re_i_17q = re.compile(re_slug_before + r'i\(17\)\(q')

# inv(3) 
re_inv_3 = re.compile(r'inv\(3\)')

# inv(3)(q21q26)
re_inv_3_q21q26 = re.compile(r'inv\(3\)\(q21\.?\d{0,2}?q26')

# inv(16) 
re_inv_16 = re.compile(r'inv\(16\)')

# inv(16)(p13q22)
re_inv_16_p13q22 = re.compile(r'inv\(16\)\(p13\.?\d{0,2}?q22')

# t(3;3)
re_t_3_3 = re.compile(r't\(3;\s?3\)')

# t(6;9)
re_t_6_9 = re.compile(r't\(6;\s?9\)')

# t(6;11)
re_t_6_11 = re.compile(r't\(6;\s?11\)')

# t(8;21)
re_t_8_21 = re.compile(r't\(8;\s?21\)')

# t(9;11)
re_t_9_11 = re.compile(r't\(9;\s?11\)')

# t(9;22)
re_t_9_22 = re.compile(r't\(9;\s?22\)')

# t(15;17) 
re_t_15_17 = re.compile(r't\(15;\s?17\)')

# t(16;16)
re_t_16_16 = re.compile(r't\(16;\s?16\)')

#-----------------------------------------------------------------------------
# See knotes.txt for details and references.

