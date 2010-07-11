import os
import sys
import logging

import wsgiref.handlers
from google.appengine.ext import webapp
from google.appengine.ext.webapp import template

sys.path.insert( 0, os.path.join(
	os.path.dirname(os.path.realpath(__file__)), 'lib' ))

from genomecommons.variantanalyzer import VariantAnalyzer

class VAT(webapp.RequestHandler):
	def get(self):
		self.response.out.write( self._render() )

	def post(self):
		self.response.out.write( self._render() )

	def _render(self):
		pv = { 'pagetitle': 'Variant Analysis Tool',
			   'path': self.request.path,
			   'url': self.request.url
			   }

		varspec = self.request.get('varspec')
		if varspec is not None and varspec != '':
			pv['varspec'] = varspec
			pv['pagetitle'] = pv['pagetitle'] + ' - ' + varspec
			pv['va'] = VariantAnalyzer(varspec)
		tmpl = os.path.join( os.path.dirname(__file__),
							 'templates',
							 'vat.html' )
		return template.render( tmpl, pv )


	def _env(self):
		return map(lambda k: '<br>%s=%s\n' % (k,os.environ[k]), sorted(os.environ.keys()))
