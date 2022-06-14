#!/usr/bin/python

import subprocess
import os
import sys

GENERATE_TABLE = './backend/generate_bit_transition_cache.out'


limit = 30
k_limit = 10
if(len(sys.argv) > 1):
	limit = int(sys.argv[1])
if(len(sys.argv) > 2):
	k_limit = int(sys.argv[2])

for n in range(limit):
	for k in range(k_limit):
		if (2*n)+k > limit:
			continue
		print(' '.join([GENERATE_TABLE, str(n), str(k)]))
		subprocess.Popen([GENERATE_TABLE, str(n), str(k)]).wait()
