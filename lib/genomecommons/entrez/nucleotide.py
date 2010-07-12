import os
import re
import string
import warnings

import genomecommons.euhelpers as euhelpers

class Nucleotide(object):
	def __init__(self,gi):
		self.id = gi
		self.seqrecord = euhelpers.fetch_nuc_gi_as_seqrecord(gi)

	@property
	def url(self):
		return 'http://www.ncbi.nlm.nih.gov/nuccore/' + self.id

