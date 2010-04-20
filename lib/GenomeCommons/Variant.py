# GenomeCommons/VarSpec -- specification of a genomic variant
# This code implements the specification at http://www.hgvs.org/mutnomen/

from GCException import GCException

class GCVariantSyntaxException(GCException):
	def __init__(self): pass

class Variant(object):
	"""
	
	"""

	def __init__(self,varspec):
		self.varspec = varspec

	@property
	def accession(self):
		ac,sep,var = self.varspec.partition(':')
		print ac, sep, var


	def __validate(self): pass
	def __validate_syntax(self): pass
	def __validate_variant(self): pass


