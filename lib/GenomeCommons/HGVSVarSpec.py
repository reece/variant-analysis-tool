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


import re
import exceptions

class InvalidVariantSyntax(Exception):
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
		varspec_re = re.compile('^([^:]+):([gp])\.(.+)')
		m = varspec_re.match(self.varspec)
		if m == None:
			raise InvalidVariantSyntax
		self.accession,self.type = m.group(1,2)
		self.varlist = m.group(3).split(',')

	def __str__(self):
		return self.varspec
