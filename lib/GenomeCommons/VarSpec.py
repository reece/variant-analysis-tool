# GenomeCommons/VarSpec -- specification of a sequence variants
# This code implements the specification at http://www.hgvs.org/mutnomen/
# 2010-04-19 19:17 Reece Hart <reece@harts.net>

import re
import Exception

from Exception import Exception
class InvalidVariantSyntax(Exception):
	def __init__(self): pass

class VarSpec(object):
	"""
GenomeCommons.VarSpec -- specification of sequence variation
Implements the specification at http://www.hgvs.org/mutnomen/
	"""

	def __init__(self,varspec):
		self.varspec = varspec
		self.accession = None
		self.type = None
		self.variantlist = None

		# TODO: The following regexp matches a subset of the HGVS variant spec
		
		varspec_re = re.compile('^([^:]+):([cgmpr])\.(.+)')
		m = varspec_re.match(self.varspec)
		if m == None:
			raise InvalidVariantSyntax
		self.accession,self.type = m.group(1,2)
		self.varlist = m.group(3).split(',')

	def __str__(self):
		return self.varspec
