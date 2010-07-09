# CoordinateMapper -- map between genomic, cds, and protein coordinates
# AUTHOR: Reece Hart <reecehart@gmail.com>
# LICENSE: BioPython

# Examples:
# AB026906.1:g.7872G>T
# AB026906.1:c.274G>T
# BA...:p.Asp92Tyr
# 
# Each refers to congruent variants in different sequence types and
# coordinates. A coordinate mapper is needed to translate between at least
# these three coordinate frames.  The mapper should deal with specialized
# syntax for splicing and UTR (e.g., 88+1, 89-2, -14, *46). In addition,
# care should be taken to ensure consistent 0- or 1-based numbering (0
# internally, as with Python/BioPython and Perl/BioPerl).
# 
# genomic  ------00000000000-----1111111----------22222222222*222-----
#          0     s0        e0    s1    e1         s2            e2
#                \         \     |     |          /             /
#                 +--+      +--+ |     | +-------+     +-------+          
#                     \         \|     |/             /
# cds                 00000000000111111122222222222*222
#                     0         c0     c1             c2
# codons:             aaabbbcccdddeeefffggghhhiiijj*kkk
# protein:              A  B  C  D  E  F  G  H  I  J  K
#                       0  p1 p2 ...                  pn
# 
# N.B. coordinates above and in code are "Python counting" (i.e., 0 based,
# right open interval).
#
# TODO:  
# * implement extended syntax for cds coords (see HGVS spec). Briefly,
# these enable representions in non-coding regions, such as 5' UTR (eg,
# -14), 3' UTR (*46), or introns (88+1 or 89-2).
# * modify g2c and c2g to use the extended sytax
# * use AbstractLocation & friends in lieu of integer locations, or at
# least accept them for initialization?  Perhaps add specialized class for
# extended syntax. Also consider __str__ implementation for ranges.
# * consider making all conversions as list-to-list. This would allow
# users to specify ranges or enumerations in any coordinate, and compute
# the resulting range/enumeration in another.  Single positions become
# degenerate cases of this structure.
# * need lots of error handling for bogus coords (e.g., cpos outside cds)

class CoordinateMapper(object):
    def __init__(self,selist=None,seqrecord=None):
        if seqrecord is not None:
            self.exons = self._extract_cds_from_seqrecord(seqrecord)
        else:
            self.exons = selist

    def _extract_cds_from_seqrecord(self,sr):
        cdsf = [ f for f in sr.features if f.type == 'CDS' ][0]
        return [ (sf.location.nofuzzy_start,
				  sf.location.nofuzzy_end)
                 for sf in cdsf.sub_features ]

    def genome_to_cds(self,gpos):
        d = 0
        for s,e in self.exons:
            l = e - s
            if gpos < s:
                return None
            if gpos <= e:
                return d+gpos-s
            d += l
        return None

    def cds_to_genome(self,cpos):
        d = 0
        for s,e in self.exons:
            l = e - s
            if cpos < d+l:
                return s+cpos-d
            d += l
        return None

    def cds_to_protein(self,cpos):
        return int(cpos/3)

    def protein_to_cds(self,ppos):
        return ppos*3,ppos*3+2          # inclusive. use python counting?



if __name__ == '__main__':
    # The following exons are from AB026906.1.
    # test case: g.7872 -> c.274 -> p.92
    # From http://www.mutalyzer.nl/1.0.4/
    exons = [ (5808,5860), (6757,6874), (7767,7912), (13709,13785) ]
    cm = CoordinateMapper(exons)
    for g1 in (7870,7871,7872,7873,7874):
        c1 = cm.genome_to_cds(g1)
        p1 = cm.cds_to_protein(c1)
        c2 = cm.protein_to_cds(p1)
        g2 = cm.cds_to_genome(c2[0])
        print('g.%s -> c.%s -> p.%s -> c.%s -> g.%s'
              % (g1+1,c1+1,p1+1,c2,g2+1))
