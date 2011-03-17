import os, re, string, warnings

from Bio import Entrez

import euhelpers as euhelpers
from hgvsvarspec import HGVSVarSpec
from coordinatemapper import CoordinateMapper
from entrez.gene import Gene
from entrez.nucleotide import Nucleotide


Entrez.email = 'reece@berkeley.edu'
Entrez.tool = __file__
os.linesep = '\n'


class VariantAnalyzer(object):
	def __init__(self,varspec):
		self._derived_varspecs = 0
		self.orig_varspec = HGVSVarSpec(varspec)
		self.g_varspec = None
		self.c_varspec = None
		self.p_varspec = None
		if self.orig_varspec.type == 'c':
			c_varspec = self.orig_varspec
		elif self.orig_varspec.type == 'g':
			g_varspec = self.orig_varspec
		elif self.orig_varspec.type == 'p':
			p_varspec = self.orig_varspec
		else:
			raise ValueError('Variant type %s is not implemented'
							 % self.orig_varspec.type)

	@property
	def gene(self):
		try:
			return self._gene
		except AttributeError:
			self._gene = Gene(self.gene_id)
			return self._gene

	@property
	def omims(self):
		return euhelpers.fetch_omim_records(
			euhelpers.link_gene_id_to_omim_ids(self.gene_id))

	@property
	def pubmeds(self):
		return euhelpers.fetch_pubmed_records(
			euhelpers.link_gene_id_to_pubmed_ids(self.gene_id))

	@property
	def variants(self):
		return euhelpers.fetch_snp_records(
			euhelpers.link_gene_id_to_snp_ids(self.gene_id))



	@property
	def gene_id(self):
		return euhelpers.link_nuc_gi_to_gene_id(self.nuc_gi)

	@property
	def nuc(self):
		try:
			return self._nuc
		except AttributeError:
			self._nuc = Nucleotide(self.nuc_gi)
			return self._nuc

	@property
	def nuc_ac(self):
		if self.orig_varspec.type in ('c','g'):
			return self.orig_varspec.accession
		return None

	@property
	def nuc_gi(self):
		return euhelpers.link_ac_to_nuc_gi(self.nuc_ac)


	def links(self):
		#http://www.ncbi.nlm.nih.gov/gene/6392
		#http://www.ncbi.nlm.nih.gov/sites/varvu?gene=6392
		#http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?locusId=6392
		pass


	def derive_varspecs(self):
		if self._derived_varspecs != 0:
			return
		self._derived_varspecs = 1

		ovs = self.orig_varspec
		m = re.match(r'(\d+)(\w>\w)',ovs.varlist[0])
		if m is None:
			raise SyntaxError
		pos = int(m.group(1))
		mut = m.group(2)

		cm = CoordinateMapper(seqrecord=self.nuc.seqrecord)
		if ovs.type == 'p':
			# no attempt made to infer cds or genomic varspec
			# c, g varspecs are undefined
			return
		elif ovs.type == 'g':
			var0 = ovs.var_i(0)
			vartxt = '%s:c.%s%s' % (ovs.accession,
				cm.genome_to_cds(var0['pos']-1)+1, var0['mut'])
			self.c_varspec = HGVSVarSpec(vartxt)
		elif ovs.type == 'c':
			var0 = ovs.var_i(0)
			print var0
			vartxt = '%s:g.%s%s' % (ovs.accession,
				cm.cds_to_genome(var0['pos']-1)+1, var0['mut'])
			self.g_varspec = HGVSVarSpec(vartxt)

		var0 = self.c_varspec.var_i(0)
		# FIXME: the following are placeholders for the real derived
		# sequence
		pro_ac = 'proid'
		aa_wt = 'Aaa'
		aa_var = 'Bbb'
		vartxt = '%s:p.%s%s%s' % (pro_ac,
					aa_wt, cm.cds_to_protein(var0['pos']-1)+1, aa_var)
		self.p_varspec = HGVSVarSpec(vartxt)



if __name__ == '__main__':
	from pprint import pprint

	for vstext in [
		'AB026906.1:g.7872G>T',
		]:
		print('* %s' % vstext)
		va = VariantAnalyzer(vstext)
		print 'ac:', va.nuc_ac
		print 'nuc_gi:', va.nuc_gi


# snp urls:
# http://www.ncbi.nlm.nih.gov/projects/SNP/SNPeutils.htm
# http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&id=3000&report=Brief
# 
