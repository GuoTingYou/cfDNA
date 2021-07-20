from ..base import ReadFileBase

class PassLocSelector:

    def __init__(self, vcfFile, pass_only=True):
        self.vcfFile = vcfFile
        self.pass_only = pass_only

    def __iter__(self):
        for variant in self.vcfFile:
            if variant.is_empty:
                break
            if self.pass_only:
                if variant.is_pass:
                    yield variant
            else:
                yield variant

ReadFileBase.register(PassLocSelector)
