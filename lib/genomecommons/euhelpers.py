# euhelpers -- Entrez Utilities helper functions
# these routines facilitate access to Entrez
# Results are memoized for speed

# http://www.ncbi.nlm.nih.gov/entrez/query/static/entrezlinks.html
# http://eutils.ncbi.nlm.nih.gov/corehtml/query/DTD/eLink_090910.dtd

import os
import sys
import shelve
import string

from Bio import Entrez
from Bio import SeqIO

from genomecommons.exceptions import *

# FIXME: make memoize do the right thing depending
# on environment (delegate entirely to memcache?)
if os.environ.has_key('SCRIPT_NAME'):
	from memoize_gae import memoize
else:
	from memoize import memoize


Entrez.email = 'reece@berkeley.edu'
Entrez.tool = __file__
os.linesep = '\n'

def ac_to_nuc_gi(ac):
	gis = ac_to_nuc_gis(ac)
	if len(gis) == 0:
		raise GCError(
			'%s: no such nucleotide record' % (ac))
	elif len(gis) > 1:
		raise GCError(
			'more than one record returned for term %s' % (ac))
	return gis[0]

def ac_to_nuc_gis(ac):
	h = Entrez.esearch(db='nucleotide',term=ac)
	r = Entrez.read(h)
	h.close()
	return r['IdList']

def fetch_nuc_gi_as_seqrecord(gen_gi):
	h = Entrez.efetch(db='nucleotide', id=gen_gi, rettype='gb')
	return SeqIO.parse(h,'genbank').next()

def fetch_gene_id_as_dict(gene_id):
	return Entrez.read(Entrez.efetch(
		db='gene',id=gene_id,retmode='xml'))[0]

def link_nuc_gi_to_gene_id(gen_gi):
	# TODO: is there ever more than one gene_id per gi?
	r = Entrez.read(Entrez.elink(
		dbfrom='nucleotide',db='gene',id=gen_gi,retmode='xml'))
	return r[0]['LinkSetDb'][0]['Link'][0]['Id']



if __name__ == '__main__':
	ac = 'AB026906.1'
	print 'ac:', ac
	print 'ac_to_nuc_gis:', ac_to_nuc_gis(ac)
	print 'ac_to_nuc_gi:', ac_to_nuc_gi(ac)
