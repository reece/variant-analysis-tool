import os
import sys
import logging

import wsgiref.handlers
from google.appengine.ext import webapp
from google.appengine.ext.webapp import template

sys.path.insert( 0, os.path.join(
	os.path.dirname(os.path.realpath(__file__)), 'lib' ))

from genomecommons.variantanalyzer import VariantAnalyzer

template_root = os.path.join( os.path.dirname(__file__), 'templates' )

class Redirect(webapp.RequestHandler):
	def get(self):
		self.redirect("/vat")
	def post(self):
		self.redirect("/vat")


class VAT(webapp.RequestHandler):
	template_dir = template_root

	def get(self):
		self.response.out.write( self._render() )

	def post(self):
		self.response.out.write( self._render() )

	def _env(self):
		return map(lambda k: '<br>%s=%s\n'
				   % (k,os.environ[k]), sorted(os.environ.keys()))

	def template_path(self):
		return os.path.join( self.template_dir,
							 os.path.basename(self.request.path) + '.html' )

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
		return template.render( self.template_path(), pv )


class TabHandler(VAT):
	template_dir = os.path.join( template_root, 'tab' )

	def _render(self):
		varspec = self.request.get('varspec')
		if varspec is not None and varspec != '':
			pv = { 'va': VariantAnalyzer(varspec),
				   'env' : self._env() }
			return template.render( self.template_path(), pv )
		else:
			return 'varspec query argument not provided'
		

