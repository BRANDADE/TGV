class DisplayRow:
    def __init__(self):
        self.locus = None
        self.depth = None
        self.size = None
        self.motifs = None

        self.genotype = None
        self.classification = None
        self.rep1 = None
        self.rep2 = None
        self.seg1 = None
        self.seg2 = None



class DisplayDetails:
    def __init__(self):
        self.locus = None
        self.depth = None
        self.size = None
        self.motifs = None
        self.purity = None
        self.methylation = None

        self.classification = None
        
        self.sequence1 = None
        self.sequence2 = None
        self.interruptions1 = None
        self.interruptions2 = None
        self.motifs_use1 = None
        self.motifs_use2 = None
        self.seg1 = None
        self.seg2 = None
        self.rep1 = None
        self.rep2 = None


class DisplayExport:
    def __init__(self):
        self.locus = None
        self.depth1 = None
        self.depth2 = None
        self.motifs = None

        self.genotype = None
        self.classification = None



class DisplayHtml:
    def __init__(self):
        self.locus = None
        self.depth = None
        self.motifs = None
        self.genotype = None
        self.classification = None
