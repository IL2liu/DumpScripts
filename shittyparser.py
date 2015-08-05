#! -*- coding:utf-8- -*-

import re
import sys
from gzopen import gzopen

with gzopen(sys.argv[1]) as f:
   for line in f:
      shit,score = line.rstrip('\n').split('\t')
      if int(score) < 11: continue
      pair = re.sub("[')(]", '', shit).replace(' ', '_').split(',_')
      sys.stdout.write('%s (u) %s = %s\n' % (pair[0], pair[1], score))
      #sys.stdout.write('%s u %s\n' % tuple(pair))
