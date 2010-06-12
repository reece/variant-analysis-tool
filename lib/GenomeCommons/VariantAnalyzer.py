import os
import re
import shelve

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
		self.vs.validate()

	def print_summary(self):
		vs = self.vs
		print('%s = %s / %s / %s' % (vs, vs.accession, vs.type, vs.varlist))

		gi = ac_to_gi(vs.accession)
		# TODO: no such!

		r = gi_as_seqrecord(gi)

		print('id:%s, gi:%s' % (r.id, gi))
		print('description:%s' % (r.description))
		print('features:%d' % (len(r.features)))
		#for f in r.features:
		#	print('####################')
		#	print(f)

		pos = int( re.search('^(\d+)',vs.varlist[0]).group(0) )
		seq = r.seq
		try:
			print vs.varlist[0], pos, seq[pos-1]
		except:
			print r
