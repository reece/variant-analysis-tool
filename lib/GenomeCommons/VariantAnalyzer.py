import os
import re
import shelve
import warnings

import Bio.Entrez
import Bio.SeqIO

import GenomeCommons.HGVSVarSpec
from GenomeCommons.utils import *


Bio.Entrez.email = 'reece@berkeley.edu'
Bio.Entrez.tool = __file__
#os.linesep = '\n'


class VariantAnalyzer(object):
	def __init__(self,vartxt):
		self.vs = GenomeCommons.HGVSVarSpec.HGVSVarSpec(vartxt)

	def print_summary(self):
		vs = self.vs
		vs.validate()

		print('%s = %s / %s / %s' % (vs, vs.accession, vs.type, vs.varlist))

		gi = ac_to_gi(vs.accession)
		r = gi_as_seqrecord(gi)
		print('id:%s, gi:%s' % (r.id, gi))
		print('description:%s' % (r.description))
		print('features:%d (%s)'
			  % (len(r.features), ','.join([f.type for f in r.features])))


		seq = r.seq
		pos = int( re.search('^(\d+)',vs.varlist[0]).group(0) )
		print vs.varlist[0], pos, seq[pos-1]
		# TODO assert that subst na == seq[pos-1] 
		
		print '- '*20
		cdsf = [ f for f in r.features if f.type == 'CDS' ][0]
		print 'cds feature:', cdsf
		print '* %s %s (%s)' % (f.type,f.location,seq[0:9])
		for sf in cdsf.sub_features:
			seq = 'NA'
			try:
				seq = sf.extract(r).seq
			except:
				pass
			print '** %s %s (%s)' % (sf.type,sf.location,seq[0:9])
		for gi in [5809,5810,5811,
				   5812,5813,5814,
				   7872
				   ]:
			ci = gidx_to_cidx(cdsf,gi-1)
			if ci is None:
				ci = 'NC'
				pi = 'NC'
			else:
				pi = cidx_to_pidx(ci)
			print('%d -> %s -> %s' % (gi,ci+1,pi+1))
