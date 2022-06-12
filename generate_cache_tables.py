#!/usr/bin/python

import subprocess
import os
import sys

GENERATE_TABLE = './backend/generate_bit_transition_cache.out'
path = './backend/transition_counts'
header_template = 'count_%d_%d.h'
source_template = 'count_%d_%d.cc'


# case_template = 'if(n == %d and k == %d){\n\t\treturn get_transition_count_%d_%d(trans_num, rec_num);\n\t}'
case_template = '\t\t\tcase %d:\n\t\t\t\treturn get_transition_count_%d_%d(trans_num, rec_num);\n'
switch_k_template = '\t\tcase %d:\n\t\tswitch(k){\n'

includes = ['#pragma once']

cc_includes = []

func = 'inline size_t get_transition_count_cache(size_t n, size_t trans_num, size_t k, size_t rec_num){\n\tswitch(n){\n'

limit = 30
k_limit = 10
if(len(sys.argv) > 1):
	limit = int(sys.argv[1])
if(len(sys.argv) > 2):
	k_limit = int(sys.argv[2])

for n in range(limit):
	func += switch_k_template % n
	for k in range(k_limit):
		if (2*n)+k > limit:
			continue
		print(' '.join([GENERATE_TABLE, str(n), str(k), os.path.join(path, source_template % (n, k)), 
			os.path.join(path, header_template % (n, k)), header_template % (n, k)]))
		subprocess.Popen([GENERATE_TABLE, str(n), str(k), os.path.join(path, source_template % (n, k)), 
			os.path.join(path, header_template % (n, k)), header_template % (n, k)]).wait()
		includes.append('#include "%s"' % (header_template % (n,k)))
		cc_includes.append('#include "%s"' % (source_template % (n,k)))

		func += case_template % (k, n, k)
	func += '\t\t}\n\t\tbreak;\n'

func += '\t}\n\tfprintf(stderr, "Unimplemented n,k tuple %lu, %lu", n, k); exit(4);\n}\n'

# case_code = ' else '.join(cases)
# func = 'inline size_t get_transition_count_cache(size_t n, size_t trans_num, size_t k, size_t rec_num){\n\t' + case_code +\
# 		' else {\n\t\tfprintf(stderr, "Unimplemented n,k tuple %lu, %lu", n, k);}\n}\n'

with open(os.path.join(path, 'cached_transition_count.h'), 'w') as f:
	f.write('\n'.join(includes))
	f.write('\n')
	f.write(func)

with open(os.path.join(path, 'cached_transition_count.cc'), 'w') as f:
	f.write('\n'.join(cc_includes))
	f.write('\n')

