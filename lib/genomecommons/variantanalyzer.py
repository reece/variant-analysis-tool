import os
import re
import warnings

from Bio import Entrez

from genomecommons.utils import *
from genomecommons.hgvsvarspec import HGVSVarSpec
from coordinatemapper import CoordinateMapper

from pprint import pprint

Bio.Entrez.email = 'reece@berkeley.edu'
Bio.Entrez.tool = __file__
os.linesep = '\n'

class VariantAnalyzer(object):
	def __init__(self,varspec):
		self.origvs = HGVSVarSpec(varspec)
		self.vs = {}
		self.vs[self.origvs.type] = self.origvs	# type in {c,g,p}


	@property
	def gen_gi(self):
		try:
			return ac_to_gi(self.vs['g'].accession)
		except KeyError:
			return None
		else:
			pass

	@property
	def gen_id(self):
		try:
			return self.vs['g'].accession
		except KeyError:
			return None
		else:
			pass

	@property
	def gen_seqrecord(self):
		try:
			return self._gen_seqrecord
		except AttributeError:
			h = Bio.Entrez.efetch(db='nucleotide', id=gi, rettype='gb')
			self._gen_seqrecord = Bio.SeqIO.parse(h,'genbank').next()
			return self._gen_seqrecord

	@property
	def gene(self):
		return self.gene_record['Entrezgene_gene']['Gene-ref']['Gene-ref_locus']

	@property
	def gene_description(self):
		return self.gene_record['Entrezgene_gene']['Gene-ref']['Gene-ref_desc']

	@property
	def gene_id(self):
		r = Entrez.read(Entrez.elink(
			dbfrom='nucleotide',db='gene',id=self.gen_gi,retmode='xml'))
		# TODO: is there ever more than one?
		return r[0]['LinkSetDb'][0]['Link'][0]['Id']

	@property
	def gene_record(self):
		try:
			return self._gene_record
		except AttributeError:
			self._gene_record = Entrez.read(Entrez.efetch(
				db='gene',id=self.gene_id,retmode='xml'))[0]
			return self._gene_record

	@property
	def gene_summary(self):
		return self.gene_record['Entrezgene_summary']
	
	@property
	def gene_synonyms(self):
		return self.gene_record['Entrezgene_gene']['Gene-ref']['Gene-ref_syn']


	@property
	def links(self):
		assert(0)
		# the following is one source of links
		# finish this function when other sources are in hand
		return self.gene_record['Entrezgene_gene']['Gene-ref']['Gene-ref_db']

	@property
	def locus(self):
		return self.gene_record['Entrezgene_gene']['Gene-ref']['Gene-ref_maploc']

	@property
	def species(self):
		r = self.gene_record['Entrezgene_source']['BioSource']['BioSource_org']['Org-ref']
		return '%s (%s)' % (r['Org-ref_taxname'], r['Org-ref_common'])


#		cm = CoordinateMapper(seqrecord=r)
#
#		if self.vs.type == 'g':
#			gpos = pos1 - 1
#			cpos = cm.genome_to_cds(gpos)
#		elif self.vs.type == 'c':
#			cpos = pos1 - 1
#			gpos = cm.cds_to_genome(cpos)
#		ppos = cm.cds_to_protein(cpos)



