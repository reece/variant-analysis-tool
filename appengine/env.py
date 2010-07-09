#!/usr/bin/python

import os

print('Content-type: text/plain')
print

for k in sorted(os.environ.keys()):
	print('%s=%s' % (k, os.environ[k]))
