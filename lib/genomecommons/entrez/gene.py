import os
import re
import string
import warnings

import genomecommons.euhelpers as euhelpers

class Gene(object):
	def __init__(self,gene_id):
		self.id = gene_id
		self.gene_dict = euhelpers.fetch_gene_records(self.id)[0]

	@property
	def description(self):
		return self.gene_dict['Entrezgene_gene']['Gene-ref']['Gene-ref_desc']

	@property
	def name(self):
		return self.gene_dict['Entrezgene_gene']['Gene-ref']['Gene-ref_locus']

	@property
	def maploc(self):
		return self.gene_dict['Entrezgene_gene']['Gene-ref']['Gene-ref_maploc']

	@property
	def species(self):
		r = self.gene_dict['Entrezgene_source']['BioSource']['BioSource_org']['Org-ref']
		return '%s (%s)' % (r['Org-ref_taxname'], r['Org-ref_common'])

	@property
	def summary(self):
		return string.replace(self.gene_dict['Entrezgene_summary'],
							  '[provided by RefSeq]', '')
	@property
	def synonyms(self):
		return self.gene_dict['Entrezgene_gene']['Gene-ref']['Gene-ref_syn']

	@property
	def url(self):
		return 'http://www.ncbi.nlm.nih.gov/gene/' + str(self.id)

	@property
	def links(self):
		assert(0)
		# the following is one source of links
		# finish this function when other sources are in hand
		return self.gene_dict['Entrezgene_gene']['Gene-ref']['Gene-ref_db']



if __name__ == '__main__':
	from pprint import pprint

	g = Gene(6392)
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
gene: %(name)s (%(id)d; %(maploc)s)  
description %(desc)s
species: %(species)s
url: %(url)s
summary: %(summary)s
synonyms: %(synonyms)s
''' % data
