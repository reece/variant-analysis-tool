import os,sys
import logging

import wsgiref.handlers
from google.appengine.ext import webapp
from google.appengine.ext.webapp import template

sys.path.insert( 0, os.path.join(
	os.path.dirname(os.path.realpath(__file__)), 'lib' ))

import genomecommons.variantanalyzer

class VAT(webapp.RequestHandler):
	def get(self):
		self.response.out.write( self._render() )

	def post(self):
		self.response.out.write( self._render() )

	def _render(self):
		pv = { 'pagetitle': 'Variant Analysis Tool',
			   'path': self.request.path,
			   'sys': sys,
			   }

		vs =  self.request.get('varspec')
		if vs is not None and vs != '':
			pv['varspec'] = vs
			pv['pagetitle'] = pv['pagetitle'] + ' - ' + vs
			va = genomecommons.variantanalyzer.VariantAnalyzer(vs)

		tmpl = os.path.join( os.path.dirname(__file__),
							 'templates/vat.html' )
		return template.render( tmpl, pv )
