import sys
from collections import Counter

c = Counter(sys.stdin)
for x in sorted(c):
  print("%s\t%s" %(x.strip(), c[x]))
