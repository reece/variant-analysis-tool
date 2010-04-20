import os,sys
import logging
import wsgiref.handlers
from google.appengine.ext import webapp
from google.appengine.ext.webapp import template
from Bio import Entrez
from Bio import SeqIO

Entrez.email = 'reece@berkeley.edu'
Entrez.tool = 'Genome Commons > Variant Analysis Tool'
os.linesep = '\n'

class VAT(webapp.RequestHandler):
    def _render(self,vals):
        page_vars = { 'pagetitle': 'Variant Analysis Tool',
                      'path': self.request.path,
                      'sys': sys,
                      'acc': self.request.get('acc') }
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
        if 0: # remove when done testing
            acc = self.request.get('acc')
            vs = self.request.get('varspec')
            handle = Entrez.esearch(db='nucleotide', term=acc)
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

# fetches raw text
#        record = handle.read()
#        self.response.out.write( self._render( {
#                    'record': record
#                    } ) )

        gb = SeqIO.parse(handle,'gb').next()
        self.response.out.write( self._render( {
                    'gi': gi,
                    'id': gb.id,
                    'name': gb.name,
                    'description': gb.description
                    } ) )
        
        handle.close()
        return


# From PC:
#from Bio import Entrez
#from Bio import SeqIO
#handle = Entrez.efetch(db='protein', id='114391',rettype='gp')
#record = SeqIO.read(handle, 'genbank')
