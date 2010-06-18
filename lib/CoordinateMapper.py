# CoordinateMapper -- map between genomic, cds, and protein coordinates
# AUTHOR: Reece Hart <reecehart@gmail.com>
# LICENSE: BioPython

# Examples:
# AB026906.1:g.7872G>T
# AB026906.1:c.274G>T
# BA...:p.Asp92Tyr
# 
# All refer to congruent variants. A coordinate mapper is needed to
# translate between at least these three coordinate frames.  The mapper
# should deal with specialized syntax for splicing and UTR (e.g., 88+1,
# 89-2, -14, *46). In addition, care should be taken to ensure consistent 0-
# or 1-based numbering (0 internally, as with Python/BioPython and
# Perl/BioPerl).
# 
# g   -----------00000000000-----1111111----------22222222222*222-----
#                s0        e0    s1    e1         s2            e2
#                \         \     |     |          /             /
#                 +--+      +--+ |     | +-------+     +-------+          
#                     \         \|     |/             /
# c                   00000000000111111122222222222*222
#                               c0     c1             c2
#                     aaabbbcccdddeeefffggghhhiiijj*kkk
# p                     A  B  C  D  E  F  G  H  I  J  K
#                       p0 p1 p2 ...                  pn
# 
# 
# TODO:  
# * g2c returns index + extended syntax
# * c2g accepts index + extended syntax


class CoordinateMapper(object):
	def __init__(self,selist=None,seqrecord=None):
		if seqrecord is not None:
			self.exons = self.__extractCDSFromSeqRecord(seqrecord)
		else:
			self.exons = selist

	def __extractCDSFromSeqRecord(self,sr):
		cdsf = [ f for f in sr.features if f.type == 'CDS' ][0]
		return [ (sf.location.start.position,sf.location.end.position)
				 for sf in cdsf.sub_features ]

	def g2c(self,gpos):
		d = 0
		for s,e in self.exons:
			l = e - s
			if gpos < s:
				return None
			if gpos <= e:
				return d + gpos - s
			d += l
		return None

	def c2g(self,cpos):
		d = 0
		for s,e in self.exons:
			l = e - s
			if cpos < d+l:
				return s + cpos - d
			d += l
		return None

	def c2p(self,cpos):
		return int(cpos/3)

	def p2c(self,ppos):
		return ppos*3,ppos*3+2


if __name__ == '__main__':
	# The following exons are from AB026906.1.
	# test case: g.7872 -> c.274 -> p.92
	# N.B. These are python counting coordinates (0-based)
	exons = [ (5808,5860), (6757,6874), (7767,7912), (13709,13785) ]
	
	cm = CoordinateMapper(exons)
	for g1 in (7870,7871,7872,7873,7874):
		c1 = cm.g2c(g1)
		p1 = cm.c2p(c1)
		c2 = cm.p2c(p1)[0]
		g2 = cm.c2g(c2)
		print g1,c1,p1,' | ',c2,g2
