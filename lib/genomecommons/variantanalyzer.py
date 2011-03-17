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
		self.origvs = HGVSVarSpec(varspec)
		self.vs = { 'c': None, 'g': None, 'p': None }
		self.vs[self.origvs.type] = self.origvs	# type in {c,g,p}
		self._derived_varspecs = 0

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
		if self.origvs.type in ('c','g'):
			return self.origvs.accession
		return None

	@property
	def nuc_gi(self):
		return euhelpers.link_ac_to_nuc_gi(self.nuc_ac)


	def links(self):
		#http://www.ncbi.nlm.nih.gov/gene/6392
		#http://www.ncbi.nlm.nih.gov/sites/varvu?gene=6392
		#http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?locusId=6392
		pass

	@property
	def g_varspec(self):
		self.derive_varspecs()
		return self.vs['g']

	@property
	def c_varspec(self):
		self.derive_varspecs()
		return self.vs['c']

	@property
	def p_varspec(self):
		self.derive_varspecs()
		return self.vs['p']


	def derive_varspecs(self):
		if self._derived_varspecs != 0:
			return
		self._derived_varspecs = 1

		ovs = self.origvs
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
			self.vs['c'] = HGVSVarSpec(vartxt)
		elif ovs.type == 'c':
			var0 = ovs.var_i(0)
			print var0
			vartxt = '%s:g.%s%s' % (ovs.accession,
				cm.cds_to_genome(var0['pos']-1)+1, var0['mut'])
			self.vs['g'] = HGVSVarSpec(vartxt)

		var0 = self.vs['c'].var_i(0)
		pro_ac = 'proid'
		aa_wt = 'Aaa'
		aa_var = 'Bbb'
		vartxt = '%s:p.%s%s%s' % (pro_ac,
					aa_wt, cm.cds_to_protein(var0['pos']-1)+1, aa_var)
		self.vs['p'] = HGVSVarSpec(vartxt)



if __name__ == '__main__':
	from pprint import pprint

	vstext = 'AB026906.1:g.7872G>T'
	va = VariantAnalyzer(vstext)
	print 'nuc_ac:', va.nuc_ac
	print 'nuc_gi:', va.nuc_gi

	pprint(va.pubmeds)
