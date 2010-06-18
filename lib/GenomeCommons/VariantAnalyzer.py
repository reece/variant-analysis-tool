import os
import re
import shelve
import warnings

import Bio.Entrez
import Bio.SeqIO

import GenomeCommons.HGVSVarSpec
from GenomeCommons.utils import *

from CoordinateMapper import CoordinateMapper


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

		cm = CoordinateMapper(seqrecord=r)
		pos1 = int( re.search('^(\d+)',vs.varlist[0]).group(0) )
		if vs.type == 'g':
			gpos = pos1 - 1
			cpos = cm.g2c(gpos)
		elif vs.type == 'c':
			cpos = pos1 - 1
			gpos = cm.c2g(cpos)
		ppos = cm.c2p(cpos)
		print('%d -> %s -> %s' % (gpos+1,cpos+1,ppos+1))

		#seq = r.seq
		#print vs.varlist[0], pos, seq[pos-1]
		# TODO assert that subst na == seq[pos-1] 
		
