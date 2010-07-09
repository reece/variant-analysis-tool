import os
import re
import shelve
import warnings

import Bio.Entrez
import Bio.SeqIO

import genomecommons.hgvsvarspec
from genomecommons.utils import *

from coordinatemapper import CoordinateMapper


Bio.Entrez.email = 'reece@berkeley.edu'
Bio.Entrez.tool = __file__
#os.linesep = '\n'


class VariantAnalyzer(object):
	def __init__(self,vartxt):
		self.vs = genomecommons.hgvsvarspec.HGVSVarSpec(vartxt)

	def gen_gi(self):
		return ac_to_gi(self.vs.accession)
		
	def summary(self):
		gen_gi = self.gen_gi()
		r = gi_as_seqrecord(gen_gi)
		cm = CoordinateMapper(seqrecord=r)

		pos1 = int( re.search('^(\d+)',self.vs.varlist[0]).group(0) )
		if self.vs.type == 'g':
			gpos = pos1 - 1
			cpos = cm.genome_to_cds(gpos)
		elif self.vs.type == 'c':
			cpos = pos1 - 1
			gpos = cm.cds_to_genome(cpos)
		ppos = cm.cds_to_protein(cpos)

		return {
			'varspec': self.vs.varspec,
			'vs': { 'ac': self.vs.accession },
			'gen': { 'gi': self.gen_gi(),
					 'description': r.description,
					 'id': r.id,
					 'pos': gpos,
				   },
			'cds': { 'pos': cpos
				   },
			'pro': { 'pos': ppos
				   }
			}
