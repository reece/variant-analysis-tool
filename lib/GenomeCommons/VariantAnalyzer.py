import os
import shelve

import Bio.Entrez
import Bio.SeqIO

import GenomeCommons.Exception as Exception
import GenomeCommons.HGVSVarSpec
from GenomeCommons.utils import ( ac_to_gi, fetch_gi_as_seq )


Bio.Entrez.email = 'reece@berkeley.edu'
Bio.Entrez.tool = __file__
#os.linesep = '\n'


class VariantAnalyzer(object):
	def __init__(self,vartxt):
		self.vs = GenomeCommons.HGVSVarSpec.HGVSVarSpec(vartxt)

	def print_summary(self):
		vs = self.vs
		print '%s = %s / %s / %s' % (vs, vs.accession, vs.type, vs.varlist)

		gi = ac_to_gi(vs.accession)
		r = fetch_gi_as_seq(gi)

		print 'id:%s, gi:%s' % (r.id, gi)
		print 'description:%s' % (r.description)
		print '#features:%d' % (len(r.features))
		for f in r.features:
			print '####################'
			print f
		
