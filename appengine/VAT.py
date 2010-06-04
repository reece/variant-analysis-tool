import os,sys
import logging
import wsgiref.handlers
from google.appengine.ext import webapp
from google.appengine.ext.webapp import template

from Bio import Entrez
from Bio import SeqIO

sys.path.insert( 0, os.path.join(
	os.path.dirname( os.path.realpath( __file__ ) ), 'lib' ))
from GenomeCommons.VarSpec import VarSpec

Entrez.email = 'reece@berkeley.edu'
Entrez.tool = 'Genome Commons > Variant Analysis Tool'
os.linesep = '\n'

class VAT(webapp.RequestHandler):
	def _render(self,vals):
		page_vars = { 'pagetitle': 'Variant Analysis Tool',
					  'path': self.request.path,
					  'sys': sys,
					  'varspec': self.request.get('varspec') }
		if vals != None:
			page_vars.update(vals)
		temp = os.path.join(
			os.path.dirname(__file__),
			'templates/VAT.html')
		return template.render(
			temp, page_vars)


	def get(self):
		self.response.out.write( self._render(None) )

	def post(self):
		vs = VarSpec( self.request.get('varspec') )
		# TODO: extend to non-nucleotide variants (should support g,c,p,r at least)
		handle = Entrez.esearch(db='nucleotide', term=vs.accession)
		record = Entrez.read(handle)
		handle.close()
		if int(record['Count']) > 1:
			self.response.out.write(
				self._render( { 'error': 'Error: '+record['Count']+' records returned.' } )
				)
			return
		gi = record['IdList'][0]
		gi = 238018044

		handle = Entrez.efetch(db='nucleotide', id=gi, rettype='gb')
		gb = SeqIO.parse(handle,'genbank').next()
		
		data = { 'gb': gb, 'gi': gi, 'vs': vs, 'description': gb.description }
		self.response.out.write( self._render( data ) )
		
		handle.close()
		return


# From PC:
#from Bio import Entrez
#from Bio import SeqIO
#handle = Entrez.efetch(db='protein', id='114391',rettype='gp')
#record = SeqIO.read(handle, 'genbank')
