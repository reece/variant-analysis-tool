# Adapted from
# http://wiki.python.org/moin/PythonDecoratorLibrary#Memoize

import logging
import os
import shelve
import sys

logger = logging.getLogger('memoize')
if 1:
	logger.setLevel(logging.DEBUG)
	ch = logging.StreamHandler()
	logger.setLevel(logging.DEBUG)
	formatter = logging.Formatter('%(message)s','%(asctime)s')
	ch.setFormatter(formatter)
	logger.addHandler(ch)
	logger.info('logger engaged, sir')

class memoize(object):
	"""Decorator that caches a function's return value each time it is called.
	If called later with the same arguments, the cached value is returned, and
	not re-evaluated.
	"""

	def __init__(self, func):
		self.func = func
		dir='/tmp/memoize-%d' % os.getuid()
		if not os.path.exists(dir):
			os.mkdir(dir,0700)			# Is there a EAFP way to do this?
		cache_fn = os.path.join(dir,'%s.cache' % func.func_name)
		self.cache = shelve.open(cache_fn)
		logger.debug('opened cache @ '+cache_fn)

	def __del__(self):
		logger.debug("closing cache")
		self.cache.close()
	  
	def __call__(self, *args):
		Ikey = str(hash(frozenset(args)))
		key = str(frozenset(args))
		if not self.cache.has_key(key):
			logger.debug('miss:'+key)
			self.cache[key] = value = self.func(*args)
		else:
			logger.debug('hit:'+key)
		return self.cache[key]
		
	def __repr__(self):
		"""Return the function's docstring."""
		return self.func.__doc__



if __name__ == "__main__":
	@memoize
	def fib(n):
		return (n > 1) and (fib(n - 1) + fib(n - 2)) or 1
	print fib(19)

