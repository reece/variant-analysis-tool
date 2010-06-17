# http://www.ncbi.nlm.nih.gov/entrez/query/static/entrezlinks.html
# http://eutils.ncbi.nlm.nih.gov/corehtml/query/DTD/eLink_090910.dtd

import os
import sys
import shelve
import string

import Bio.Entrez
import Bio.SeqIO
from GenomeCommons.Exceptions import *

from memoize import memoize

Bio.Entrez.email = 'reece@berkeley.edu'
Bio.Entrez.tool = __file__

@memoize
def ac_to_gis(ac):
	h = Bio.Entrez.esearch(db='nucleotide',term=ac)
	r = Bio.Entrez.read(h)
	h.close()
	return r['IdList']

@memoize
def ac_to_gi(ac):
	gis = ac_to_gis(ac)
	if len(gis) == 0:
		raise GCError(
			'%s: no such nucleotide record' % (ac))
	elif len(gis) > 1:
		raise GCError(
			'more than one record returned for term %s' % (ac))
	return gis[0]

@memoize
def gi_as_seqrecord(gi):
	h = Bio.Entrez.efetch(db='nucleotide', id=gi, rettype='gb')
	return Bio.SeqIO.parse(h,'genbank').next()


def get_complete_cds(record):
	"""
	record should be an instance of Bio.SeqRecord.Record
	TODO: throw exception if not a coding sequence
	"""
	for f in record.features:
		if f.type == 'CDS':
			return f.extract(record.seq)
	return None

def gidx_to_cidx(cdsf,pos):
	"""
	given a CDS parent feature (with collection of CDS subfeatures) and
	index in the genomic reference, return a CDS index (0-based).
	"""
	if cdsf.type != 'CDS':
		raise GCError('arg is not a CDS feature')
	d = 0
	for sf in cdsf.sub_features:
		s = sf.location.start.position
		e = sf.location.end.position
		l = e - s
		#print('d=%d [%d,%d] l=%d' % (d,s,e,l))
		if pos < s:
			# pos is less than the start of this exon => non-coding
			return None
		if pos <= e:
			return d + pos-s
		d += l
	return None

def cidx_to_pidx(cidx):
	return int(cidx/3)



#def gi_to_refseq_ac(gi):
	

#def elink(dbfrom,db,id):
#	r = Entrez.read(Entrez.elink(dbfrom=dbfrom,db=db,id=id))

#r = Entrez.read( Entrez.elink(dbfrom="nucleotide",db='gene',id=5295993) )
#r[0]['LinkSetDb'][0]['Link'][0]['Id']

# a remarkably complete and convoluted record:
# r= Entrez.read( Entrez.efetch(db='nucleotide',id=5295993,rettype='xml') )
