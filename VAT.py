import os
import logging
import wsgiref.handlers
from google.appengine.ext import webapp
from google.appengine.ext.webapp import template
from Bio import Entrez
from Bio import SeqIO

Entrez.email = 'reece@berkeley.edu'
Entrez.tool = 'Genome Commons > Variant Analysis Tool'

class VAT(webapp.RequestHandler):
    def _render(self,vals):
        page_vars = { 'pagetitle': 'Variant Analysis Tool',
                      'path': self.request.path,
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
        record = handle.read()
        self.response.out.write( self._render( {
                    'record': record
                    } ) )
        handle.close()

#        handle = Entrez.efetch(db='nucleotide', id=gi, rettype='xml')
#        r = Entrez.read(handle)
#        .next()
#                    'gi': gi,
#                    'name': record.name,
#                    'description': record.description,
#                    'record': record
#                    } ) )
        handle.close()

        return
