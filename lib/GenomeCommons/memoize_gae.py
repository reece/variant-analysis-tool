import logging
import sys

from google.appengine.api import memcache

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

	def __call__(self, *args):
		key = str(args)
		rv = not memcache.get(key)
		if rv is None:
			rv = self.func(*args)
			if not memcache.add(key,rv,3600):
				logging.error('Memcache set failed.')
		return rv
		

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

