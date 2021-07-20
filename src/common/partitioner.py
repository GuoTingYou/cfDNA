from pandas import read_table, DataFrame

from ..base import GenomeBase, GenomeMeta

from .gene import GeneTable

class RegionPartitioner:
    """
    Objective
    ---------
    This class is used to partition the whole genome into several smaller regions
    according to the parameter `criterion` and `requirement` when the instance of
    this class is initialized. Valid options for `criterion` is **("gene", "fixed",
    "dynamic")** and `requirement` is an INTERGER. These two parameters together
    defined how to divide a region and should be provided in the format of
    ***criterion<STR>:requirement<INT>***.

    Parameter
    ---------
    ** criterion ** : <str> Criterion for region partitioning.
    >Valid options for `criterion` is **("gene", "fixed", "dynamic")**. Refered to
    the below comments for detailed description of each option.

    ** requirement ** : <int> A region is definitive when the requirement is met
    according to the given `criterion`.
    >This parameter is only valid together with `criterion`.
    e.g. "gene:150" means 150 bp upstream transcription start site and downstream
         transcription end site is a definitive region.
    e.g. "fixed:1000000" means divided the genome every 1000000 bp.
    e.g. "dynamic:1000" means a region should ultimately cumulate 1000 loci when
         a region is defined.

    Return
    ------
    RegionPartitioner object. This object is callable, and argument for `iterator`
    should be either `ReadSNP` or `ReadMAF` object. When RegionPartitioner intance
    is called, it will return a generator of divided regions of type `_Region`.

    Example
    -------
    >>>
    partition = RegionPartitioner('gene', 1000, genetable='path/to/genetable')
    snpFile = readSNP('path/to/sample.snp.chrN.gz')
    for region in partition(snpFile, skip_index=True):
        ... # do something
    """
    valid_criterion = (
        'bed',     # <bed:STR> Requirement should be a bed file.
        'tss',     # <tss:INT> Promoter region is defined as N bp upstream and
                   # downstream of transcription start site (TSS) of a gene.
                   # In this case `genetable` should be specified.
        'gene',    # <gene:INT> Flanking gene region. The following `genetable`
                   # is default gene set which could be changed by command line.
                   # gene region is defined by the genetable's **start** and
                   # **end** fields flanked by a given INT **extention** upstream
                   # and downtream respectively. ext__start__gene_body__end__ext

        'fixed',   # <fixed:INT> Fixed region size. Divided whole genome into
                   # several regions whose size are the same as given INT.

        'dynamic', # <dynamic:INT> A region is defined when the number of loci in
                   # the region has cumulated to reach the given INT. Thus, the
                   # genome is dynamically divided.
    )
    def __init__(self, criterion, requirement, genetable=None):
        self.criterion = criterion
        self.genetable = genetable
        self._process(requirement)
        self.regions = []

    def __call__(self, iterator, skip_indel):
        self.getValue = {
            'ReadSNP': lambda locus: locus.depth_of(locus.interested_allele),
            'ReadMAF': lambda locus: locus.minor_allele_frequency,
        }[type(iterator).__name__]

        loci = self.loci_generator(iterator, skip_indel)

        if self.regions:
            return self.along_partitioner(loci)
        return {
            'bed':     self.bed_partitioner,
            'tss':     self.gene_partitioner,    # tss region.
            'gene':    self.gene_partitioner,    # gene region. <STR:INT>
            'fixed':   self.fixed_partitioner,   # fixed region size. <INT>
            'dynamic': self.dynamic_partitioner, # dynamic region size. <INT>
        }[self.criterion](loci)

    def bed_partitioner(self, loci):
        bed = self.requirement
        n_regions = bed.shape[0]
        self.regions = [r.typed for r in map(
            _Region, bed.chrom, bed.start, bed.end, bed.symbol,
            self._nested_list(n_regions), self._nested_list(n_regions)
        )]
        return self.along_partitioner(loci)

    def gene_partitioner(self, loci):
        gt = self.genetable
        n_genes = gt.shape[0]
        self.regions = list(map(
            _Region, gt.chrom, gt.start, gt.end, gt.symbol,
            self._nested_list(n_genes), self._nested_list(n_genes)
        ))
        return self.along_partitioner(loci)

    def fixed_partitioner(self, loci):
        getValue = self.getValue
        start, depths, values = None, [], []
        for locus in loci:
            end = int(locus.start / self.region_size)
            if start is None:
                start = end
            elif start != end:
                start *= self.region_size
                end   *= self.region_size
                self.regions.append(_Region(locus.chrom, start, end))
                yield _Region(locus.chrom, start, end,
                              name=None, depths=depths, values=values)
                start = None; depths.clear(); values.clear()
            depths.append(locus.depth)
            values.append(getValue(locus))
        else:
            if start is not None:
                start *= self.region_size
                end = (end + 1) * self.region_size
                self.regions.append(_Region(locus.chrom, start, end))
                yield _Region(locus.chrom, start, end,
                              name=None, depths=depths, values=values)
                start = None; depths.clear(); values.clear()

    def dynamic_partitioner(self, loci):
        getValue = self.getValue
        start, depths, values = None, [], []
        for locus in loci:
            end = locus.end
            if start is None:
                start = locus.start
            elif len(values) == self.cumulative:
                self.regions.append(_Region(locus.chrom, start, end))
                yield _Region(locus.chrom, start, end,
                              name=None, depths=depths, values=values)
                start = None; depths.clear(); values.clear()
            depths.append(locus.depth)
            values.append(getValue(locus))
        else:
            if start is not None:
                self.regions.append(_Region(locus.chrom, start, end))
                yield _Region(locus.chrom, start, end,
                              name=None, depths=depths, values=values)
                start = None; depths.clear(); values.clear()

    def along_partitioner(self, loci):
        getValue = self.getValue
        index = 0; n_regions = len(self.regions)
        region = self.regions[index]
        for locus in loci:
            while not locus.ison(region.chrom) and index < n_regions:
                region = self.regions[index]; index += 1
            if locus in region:
                region.depths.append(locus.depth)
                region.values.append(getValue(locus))
            else:
                while index < n_regions and region.end < locus.start:
                    region = self.regions[index]; index += 1
                else:
                    if locus in region:
                        region.depths.append(locus.depth)
                        region.values.append(getValue(locus))
        for region in self.regions:
            yield region
            region.depths.clear()
            region.values.clear()

    def loci_generator(self, iterator, skip_indel):
        for locus in iterator:
            if skip_indel and not locus.is_indel:
                yield locus
            else:
                yield locus

    def clear(self):
        self.regions.clear()

    def _process(self, requirement):
        if self.criterion == 'bed':
            self.requirement = read_table(requirement, header=None,
                usecols=[0, 1, 2, 3], names=['chrom', 'start', 'end', 'symbol'])

        elif self.criterion in ('tss', 'gene'):
            self.requirement = self.extention = int(requirement)
            geneTable = DataFrame(sorted(g for g in GeneTable(self.genetable)))
            geneTable.chrom = geneTable.chrom.apply(
                lambda chrom: type(chrom)(chrom.split('_')[0]))
            if self.criterion == 'tss':
                geneTable.end = geneTable.start.apply(
                    lambda start: start + self.extention)
            else:
                geneTable.end = geneTable.end.apply(
                    lambda end: end + self.extention)
            geneTable.start = geneTable.start.apply(
                lambda start: start - self.extention)
            geneTable.drop_duplicates('symbol', keep='first',
                                      inplace=True, ignore_index=True)
            self.genetable = geneTable

        elif self.criterion == 'fixed':
            self.requirement = self.region_size = int(requirement)

        elif self.criterion == 'dynamic':
            self.requirement = self.cumulative = int(requirement)

        else:
            raise ValueError(
                f'RegionPartitioner got an invalid argument {criterion} '
                f'for `criterion`, valid options are {self.valid_criterion}')

    def _nested_list(self, n):
        return [[] for _ in range(n)]

class _Region(GenomeBase, metaclass = GenomeMeta,
              expand_fields = ('name', 'depths', 'values')):
    @property
    def nloci(self):
        return len(self.values)

    @property
    def depth(self):
        return sum(self.depths)

    def __new__(cls, chrom, start, end, name=None, depths=None, values=None):
        #start = max(start, chrom.min_pos)
        #end   = min(end,   chrom.max_pos)
        name  = end - start if name in (None, '.') else name
        depths = [] if depths is None else depths
        values = [] if values is None else values
        return super().__new__(cls, chrom, start, end, name, depths, values)
