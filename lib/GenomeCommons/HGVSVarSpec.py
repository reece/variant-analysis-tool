# GenomeCommons/VarSpec -- specification of a sequence variants
# 2010-04-19 19:17 Reece Hart <reece@harts.net>

# This code current supports only a subset of the HGVS syntax at
# http://www.hgvs.org/mutnomen/.  It should aim for full implementation.

# In addition, we should consider refactoring the code to support multiple
# variant syntaxes (see below) and perhaps provide a superclass to enable
# seamless specification and translation.
# Variant standards and representations:
# - http://www.hgvs.org/mutnomen/
#   (with LRG support: http://www.lrg-sequence.org)
# - http://www.emqn.org/emqn/Mutation-Nomenclature.html
# - http://www.humgen.nl/mutalyzer/1.0.1/mutationcheck_help.html
# - BIC:
# - variant names (e.g., delta508, Factor V Leiden)

# TODO
# - more complete specification
# - subclasses generic mutations (enable other representations)
# - warning exception class (or other mechanism?)
#   - sequence type appropriate for ref seq
#   - versioned sequence
#   - c. ATG :== pos 1
#   - insertion -> dup
#   - x_yOPn => assert(y-x+1 == n)
#   - ins requires x_y

import re
import exceptions

import Bio.Alphabet

from Exceptions import *

aa_re_t = '|'.join(Bio.Alphabet.ThreeLetterProtein.letters)

class InvalidVariantSyntax(Error):
	def __init__(self): pass

class HGVSVarSpec(object):
	"""
	Bio.VarSpec.HGVSVarSpec -- specification of sequence variation
	Implements the specification at http://www.hgvs.org/mutnomen/ .
	IMPORTANT NOTES:
	1) The HGVS spec has not yet been finalized. This code will lag
	the specification. At the HVP 2010 meeting Johan den Dunnen was asked
	to finalize the specification and submit it for ISO approval.
	2) This code implements only a subset of the current specification.
	"""

	def __init__(self,varspec):
		self.varspec = varspec
		self.accession = None
		self.type = None
		self.varlist = None

		# recognize only g. and p. references with singleton
		# changes; fail on c,m,r, or compound variants
		varspec_re = re.compile('^([^:]+):([cgp])\.(.+)')
		m = varspec_re.match(self.varspec)
		if m == None:
			raise InvalidVariantSyntax
		self.accession,self.type = m.group(1,2)
		self.varlist = m.group(3).split(',')

	def __str__(self):
		return self.varspec

	def validate(self):
		if re.search('^(?:AC|AJ|AY|NM_|NP_)\d+$',self.accession):
			raise Warning(
				'%s is an unversioned reference sequence'
				% (self.accession))
		return
