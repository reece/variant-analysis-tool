# Adapted from
# http://wiki.python.org/moin/PythonDecoratorLibrary#Memoize
# See also:
# http://code.activestate.com/recipes/576642/
# http://code.activestate.com/recipes/325205-cache-decorator-in-python-24/
# http://code.activestate.com/recipes/498110-memoize-decorator-with-o1-length-limited-lru-cache/
# 


import atexit
import logging
import os
import shelve
import sys

logger = logging.getLogger('memoize')
if 0:
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
		logger.debug('opened cache for %s (%s)' % (func.func_name,cache_fn))
		atexit.register( lambda : self.cache.close() )

	def __call__(self, *args):
		#key = str(frozenset(args))
		key = str(args)
		if not self.cache.has_key(key):
			logger.debug('miss:%s(%s)' % (self.func.func_name,key))
			self.cache[key] = value = self.func(*args)
		else:
			logger.debug('hit:%s(%s)' % (self.func.func_name,key))
		return self.cache[key]
		

if __name__ == "__main__":
	@memoize
	def fib(n):
		return (n > 1) and (fib(n - 1) + fib(n - 2)) or 1
	print fib(19)

	@memoize
	def fibn(n,*argv):
		return (n > 1) and (fibn(n-1,argv) + fibn(n-2,argv)) or 1
	print fibn(19,'a','b',42,sys.argv)

	# memoization of methods is tricky in principle: do we want the same
	# result for the same args or not? Generally, not. The above memoizer
	# currently distinguishes instances via self, so different instances
	# use different cache slots. So far so good. The part I don't like
	# about this is that between process invocations, self will change: we
	# won't get a cache hit AND we'll have used file space for the orphan
	# result.
	# see http://ubuntuforums.org/showthread.php?t=1251060 for method handling
	class fibber(object):
		@memoize
		def fibber(self,n):
			return (n > 1) and (fib(n - 1) + fib(n - 2)) or 1
	i1 = fibber()
	# i1.fibber(5) fails for now because methods are memoizable with above
	i2 = fibber()
	# i2.fibber(6) fails for now

