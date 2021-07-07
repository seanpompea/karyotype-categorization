import unittest
from common import *

# Usage:
#   python -m unittest test_structural

ks = [
{'trisomy21a': '47,XX,+21' # AGT p366
,'trisomy21b': '47,XX,+21[10]/46,XX[10]'} # CA p38; mosaic.
}
]



'''
ex (from AGT)                                     numer event count     pg
47,XY,del(12)(q13),+del(12)(q13)                  1                     377


'''

class TestStructural(unittest.TestCase):

  def test_has_trisomy_at(self):
    self.assertEqual(has_trisomy_at(ks['trisomy21a'], 21), True)
    self.assertEqual(has_trisomy_at(ks['trisomy21b'], 21), True)

  def test_isupper(self):
    self.assertTrue('FOO'.isupper())
    self.assertFalse('Foo'.isupper())

  def test_split(self):
    s = 'hello world'
    self.assertEqual(s.split(), ['hello', 'world'])
    # check that s.split fails when the separator is not a string
    with self.assertRaises(TypeError):
      s.split(2)

#    def test_shouldfail(self):
#        s = 'lalala'
#        self.assertEqual(s, 'lalalooo')

if __name__ == '__main__':
  unittest.main()


