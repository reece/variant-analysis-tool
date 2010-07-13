# euhelpers -- Entrez Utilities helper functions
# these routines facilitate access to Entrez

# http://www.ncbi.nlm.nih.gov/entrez/query/static/entrezlinks.html
# http://eutils.ncbi.nlm.nih.gov/corehtml/query/DTD/eLink_090910.dtd

import os
import random
import sys
import shelve
import string
from xml.etree.ElementTree import XML

from Bio import Entrez
from Bio import SeqIO

from genomecommons.exceptions import *

rnd = random.random()
Entrez.email = 'reece+%s@berkeley.edu' % rnd
Entrez.tool = '__file__+%s' % rnd
os.linesep = '\n'

def fetch_nuc_gi_as_seqrecord(nuc_gi):
    h = Entrez.efetch(db='nucleotide',id=nuc_gi, rettype='gb')
    return SeqIO.parse(h,'genbank').next()

def fetch_gene_record(id):
    return fetch_gene_records([id])[0]
def fetch_gene_records(ids):
    return map(_remap_dict_keys,
		Entrez.read(Entrez.efetch(db='gene',id=','.join(ids),retmode='xml')))

def fetch_pubmed_record(id):
	return fetch_pubmed_records([id])[0]
def fetch_pubmed_records(ids):
    return map(_remap_dict_keys,
		Entrez.read(Entrez.efetch(db='pubmed',id=','.join(ids), retmode='xml')))

def fetch_omim_record(id):
    return fetch_omim_records([id])[0]
def fetch_omim_records(ids):
    return map(_remap_dict_keys,
		Entrez.read(Entrez.efetch(db='omim',id=','.join(ids), retmode='xml')))

def fetch_snp_record(id):
	return fetch_snp_records([id])[0]
def fetch_snp_records(ids):
    # 2010-07-12 11:57 Reece Hart <reece@harts.net> Most (all other?)
    # Entrez facilities use DTDs.  dbSNP uses XSD (with namespaces), which
    # isn't supported by Entrez.read.  Use xml.elementtree instead.
    xml = Entrez.efetch(db='snp',id=','.join(ids), retmode='xml').read()
    d = XML(xml)
    return map(_remap_dict_keys, map( _rs_elem_as_dict,
		d.findall('{http://www.ncbi.nlm.nih.gov/SNP/docsum}Rs')))


def link_ac_to_nuc_gi(ac):
    return link_ac_to_nuc_gis(ac)[0]
def link_ac_to_nuc_gis(ac):
    return Entrez.read(Entrez.esearch(db='nucleotide',term=ac))['IdList']

def link_nuc_gi_to_gene_id(nuc_gi):
    return link_nuc_gi_to_gene_ids(nuc_gi)[0]
def link_nuc_gi_to_gene_ids(nuc_gi):
    return _fetch_elinks('nucleotide','gene',nuc_gi)

def link_gene_id_to_omim_ids(gene_id):
    return _fetch_elinks('gene','omim',gene_id)

def link_gene_id_to_pubmed_ids(gene_id):
    return _fetch_elinks('gene','pubmed',gene_id)

def link_gene_id_to_snp_ids(gene_id):
    return _fetch_elinks('gene','snp',gene_id)




######################################################################
## INTERNALS

def _extract_link_ids(r):
    return [ e['Id'] for e in r[0]['LinkSetDb'][0]['Link'] ]

def _fetch_elinks(dbfrom,dbto,id):
	return _extract_link_ids( Entrez.read(Entrez.elink(
            dbfrom=dbfrom,db=dbto,id=id,retmode='xml' )) )

def _rs_elem_as_dict(rs):
    return {
        'rsId': rs.get('rsId'),
        'class': rs.get('snpClass'),
        'hgvs': map( lambda vs: vs.text, 
                     rs.findall('{http://www.ncbi.nlm.nih.gov/SNP/docsum}hgvs'))
        }

def _remap_dict_keys(d):
	"""recursively remap dict keys to [_a-z]
	WARNING: makes the terrible assumption that keys won't clash
	during remapping.
	"""
	if isinstance(d,dict):
		return dict([(_map_key(k),
					  _remap_dict_keys(d[k])) for k in d.keys()])
	elif isinstance(d,list):
		return map(_remap_dict_keys,d)
	else:
		return d

def _map_key(k):
	return k.lower().replace('-','_')



######################################################################
## Parsing dbSNP records. Two options:
# 1) with minidom
# from xml.dom import minidom
# snpxml = Entrez.efetch(db='snp', id='80338847,80338846,80338845,80338844,80338843', retmode='xml').read()
# p = minidom.parseString(snpxml)
# f = p.getElementsByTagName('Rs')
# print f[0].attributes['rsId'].nodeValue
# 2) with elementtree
# from xml.etree.ElementTree import XML
# snpxml = Entrez.efetch(db='snp', id='80338847,80338846,80338845,80338844,80338843', retmode='xml').read()
# s = XML(snpxml)
# rses = s.findall("{http://www.ncbi.nlm.nih.gov/SNP/docsum}Rs")
# print rses[0].attrib['rsId']


if __name__ == '__main__':
    from pprint import pprint

    dict1 = { u'a': 1, u'B-B': 2, u'C': { u'E': u'upper', u'F-blank-you': { u'G': u'howdy'} } }
    pprint(dict1)
    dict2 = _remap_dict_keys(dict1)
    pprint(dict2)

    ac = 'AB026906.1'
    print 'ac:', ac

    nuc_gi = link_ac_to_nuc_gi(ac)
    print 'nuc_gi:', nuc_gi

    gene_ids = link_nuc_gi_to_gene_ids(nuc_gi)
    gene_id = gene_ids[0]
    print '%d gene_ids: %s (using first)' % (len(gene_ids),','.join(gene_ids))
    
    omims = link_gene_id_to_omim_ids(gene_id)
    print '%d omims: %s' % (len(omims),','.join(omims))
    recs = fetch_omim_records(omims)
    print '%d omim recs received' % len(recs)
    pprint(recs[0])

#    pmids = link_gene_id_to_pubmed_ids(gene_id)
#    print '%d pmids: %s' % (len(pmids),','.join(pmids))
#    recs = fetch_pubmed_records(pmids)
#    print '%d pubmed recs received' % len(recs)
#
#    snpids = link_gene_id_to_snp_ids(gene_id)
#    print '%d snpids: %s' % (len(snpids),','.join(snpids))
#    recs = fetch_snp_records(snpids[:4])
#    print '%d snp recs received' % len(recs)
#    print recs
