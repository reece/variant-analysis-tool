import os
import shelve

import Bio.Entrez
import Bio.SeqIO

import GenomeCommons.HGVSVarSpec


Bio.Entrez.email = 'reece@berkeley.edu'
Bio.Entrez.tool = 'Genome Commons > Variant Analysis Tool'
#os.linesep = '\n'


class VariantAnalyzer(object):
	def __init__(self,vartxt):
		self.vs = GenomeCommons.HGVSVarSpec.HGVSVarSpec(vartxt)
		self.shelf = shelve.open('/tmp/cache')

	def __del__(self):
		self.shelf.close()
		# FIXME: can't figure out how to call super's __del__:
		# super(self.__class__,self).__del__()

	def _fetch_nucleotide(self,ac):
		skey = 'entrez/' + ac
		if not self.shelf.has_key(skey):
			handle = Bio.Entrez.esearch(db='nucleotide',term=ac)
			self.shelf[skey] = Bio.Entrez.read(handle)
			print 'fetched and wrote %s' % (skey)
			handle.close()
		else:
			print 'fetched %s from cache ' % (skey)
		return self.shelf[skey]

	def dome(self):
		print self._fetch_nucleotide(self.vs.accession)


#if int(record['Count']) > 1:
#	self.response.out.write(
#		self._render( { 'error': 'Error: '+record['Count']+' records returned.' } )
#		)
#	return
#
#gi = record['IdList'][0]

#	gi = 238018044
#
#	handle = Entrez.efetch(db='nucleotide', id=gi, rettype='gb')
#	gb = SeqIO.parse(handle,'genbank').next()
#
#	data = { 'gb': gb, 'gi': gi, 'vs': vs, 'description': gb.description }
#	self.response.out.write( self._render( data ) )
#
#	handle.close()
#	return


