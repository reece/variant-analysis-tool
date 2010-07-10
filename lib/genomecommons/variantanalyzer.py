import os
import re
import shelve
import warnings

import Bio.Entrez
import Bio.SeqIO

from genomecommons.utils import *
from genomecommons.hgvsvarspec import HGVSVarSpec
from coordinatemapper import CoordinateMapper

Bio.Entrez.email = 'reece@berkeley.edu'
Bio.Entrez.tool = __file__
#os.linesep = '\n'


## Information flow:
# Input: c., g., or p. varspec
# 3 cases for determining correspondence of g, c, p variants:
#   g->c,c->p
#   c->g,c->p
#   p->c,c->g
# After which we have g<->c<->p correspondence
# This mapping is ambiguous, an annoyance that has unknown 
# impact here.
# 
# from g:
# 
# nuc gi -> coords
#  	   -> rel snps
# 	   -> p gi
# 
# p gi -> gene
#     -> browser
#   	-> AC
#   -> strx, features, go
#   -> homology
#   -> polyphen, SIFT, others
# 
#         -> generif, pubs, OMIM, GO


class VariantAnalyzer(object):
	def __init__(self,vartxt):
		self.vs = HGVSVarSpec(vartxt)

	def gen_gi(self):
		# works only for 
		return ac_to_gi(self.vs.accession)
		
	def cds_gi(self):
		assert(0)

	def pro_gi(self):
		assert(0)

	def gene(self):
		assert(0)

	def locus(self):
		assert(0)

	def XXX(self):
		assert(0)

	def XXX(self):
		assert(0)

	def XXX(self):
		assert(0)

	def XXX(self):
		assert(0)

	def XXX(self):
		assert(0)


	#############################################################
	## summary methods
	def nuc_info(self):
		assert(0)
	def pro_info(self):
		assert(0)

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
			'varspec': self.vs,
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
