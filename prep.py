# prep.py

def snip_front(k):
  '''Remove all text from beginning up to the point where a
  two-digit (or greater) number appears -- which we assume is
  the chromosome count. 
  At this time, we knowingly remove 'mos ' and 'chi ' symbols.
  Eg mos 47,XXX[35]/45,X[15]  will become 47,XXX[35]/45,X[15]
  (example from AGT p405)
  and later will be split at the '/' into two separate items
  to analyze. .
  '''
  pass

def snip_nuc_ish(k):
  '''Remove all text from "nuc ish" onward.'''
  pass


def despace(k):
  '''Remove all spaces, except for special cases like
  "de novo" and "dir dup".'''
  pass

def snip_array(k):
  '''Remove all text from "arr" onward (that is, array/microarray
  nomenclature).'''
  pass

def snip_freetext(k):
  '''Remove freetext notes. Strategy:
  Look for the first contiguous run of 8 alphabetic characters,
  and remove all text from the beginning of the run onward.'''
  pass

def colon2semi(k):
  pass

def newline2space(k):
  pass

