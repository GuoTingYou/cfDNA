from ..base import ReadFileBase
from ..util import DoubleIterator
from .snp import SNP

class PassSnpSelector(DoubleIterator):

    def __iter__(self):
        variants = (v for v in super().__iter__() if v is not None)
        for genotype, het_variant, hom_variant in variants:
            if (het_variant.is_pass and het_variant.is_snp and
                hom_variant.is_pass and hom_variant.is_snp):
                snp = self._snp_maker(genotype, het_variant, hom_variant)
                if isinstance(snp, SNP):
                    yield snp

    def case_equal(self, variant_bg, variant_it):
        genotype = f'{variant_bg.genotype}-{variant_it.genotype}'
        if variant_bg.is_AAab(variant_it):
            return f'AA-ab:{genotype}', variant_it, variant_bg
        elif variant_bg.is_ABaa(variant_it):
            return f'AB-aa:{genotype}', variant_bg, variant_it

    def case_smaller(self, variant_bg, variant_it):
        if variant_bg.is_hetero:
            genotype = f'AB-aa:{variant_bg.genotype}-0/0'
            return genotype, variant_bg, variant_it

    def case_larger(self, variant_bg, variant_it):
        if variant_it.is_hetero:
            genotype = f'AA-ab:0/0-{variant_it.genotype}'
            return genotype, variant_it, variant_bg

    def case_iter1_unfinished(self, variant_bg):
        if variant_bg.is_hetero:
            genotype = f'AB-aa:{variant_bg.genotype}-0/0'
            return genotype, variant_bg, self.item2

    def case_iter2_unfinished(self, variant_it):
        if variant_it.is_hetero:
            genotype = 'AA-ab:0/0-{variant_it.genotype}'
            return genotype, variant_it, self.item1

    def _snp_maker(self, genotype, het_variant, hom_variant):
        if '0/0' in genotype:
            allele_A, allele_B = het_variant.gt_bases.split('/')
        else:
            allele_A = hom_variant.gt_bases.split('/')[0]
            allele_B = het_variant.gt_bases.split('/')
            try:
                allele_B.remove(allele_A)
            except ValueError:
                print((f'Warning: inconsist allele at {chrom}:{start}-{end}:'
                       f'{het_variant.gt_bases}:{hom_variant.gt_bases}'))
                return None
            else:
                allele_B = allele_B[0]

        return SNP(het_variant.chrom, het_variant.start, het_variant.end,
                   genotype, allele_A, allele_B, None, None, None, None, None)

ReadFileBase.register(PassSnpSelector)
