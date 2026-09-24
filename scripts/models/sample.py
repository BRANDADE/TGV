class Sample:
    """
    Représente un sample pour un TRID donné.
    Contient uniquement les deux allèles.
    """
    def __init__(self, name):
        self.name = name
        self.allele1 = None
        self.allele2 = None

        # GT tel qu'écrit dans le VCF (phasage éventuel, ex. "1|0")
        self.gt_raw = None

        # Classification clinique déjà calculée pour ce sample
        self.clinical_done = False

        self.result = None

    def __repr__(self):
        return f"<Sample {self.name} | GT={self.gt_raw} | A1={self.allele1} | A2={self.allele2}>"
