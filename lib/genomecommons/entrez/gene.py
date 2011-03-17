import os
import re
import string
import warnings

import genomecommons.euhelpers as euhelpers

class Gene(object):
	def __init__(self,gene_id):
		self.id = gene_id
		self.gene_dict = euhelpers.fetch_gene_record(self.id)

	@property
	def description(self):
		return self.gene_dict['entrezgene_gene']['gene_ref']['gene_ref_desc']

	@property
	def name(self):
		return self.gene_dict['entrezgene_gene']['gene_ref']['gene_ref_locus']

	@property
	def maploc(self):
		return self.gene_dict['entrezgene_gene']['gene_ref']['gene_ref_maploc']

	@property
	def species(self):
		r = self.gene_dict['entrezgene_source']['biosource']['biosource_org']['org_ref']
		return '%s (%s)' % (r['org_ref_taxname'], r['org_ref_common'])

	@property
	def summary(self):
		return string.replace(self.gene_dict['entrezgene_summary'],
							  '[provided by RefSeq]', '')
	@property
	def synonyms(self):
		return self.gene_dict['entrezgene_gene']['gene_ref']['gene_ref_syn']

	@property
	def url(self):
		return 'http://www.ncbi.nlm.nih.gov/gene/' + str(self.id)

	@property
	def links(self):
		assert(0)
		# the following is one source of links
		# finish this function when other sources are in hand
		return self.gene_dict['entrezgene_gene']['gene_ref']['gene_ref_db']



if __name__ == '__main__':
	from pprint import pprint

	g = Gene('6392')
	data = {
		'desc': g.description,
		'id': g.id,
		'maploc': g.maploc,
		'name': g.name,
		'species': g.species,
		'summary': g.summary,
		'synonyms': g.synonyms,
		'url': g.url,
		}
	print '''\
gene: %(name)s (%(id)s; %(maploc)s)  
description %(desc)s
species: %(species)s
url: %(url)s
summary: %(summary)s
synonyms: %(synonyms)s
''' % data
