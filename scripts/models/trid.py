class TRID:
    """
    Représente un locus TRGT global dans la run.
    """
    def __init__(self, trid):
        self.trid = trid

        # Position génomique
        self.chrom = None
        self.start = None
        self.end = None
        self.motifs = []
        self.motifs_rc = None

        # Config clinique (YAML)
        self.clinical = None

        # Bloc de clinical_thresholds.yaml utilisé (TRID lui-même, ou alias
        # défini dans configs/trid_aliases.yaml)
        self.clinical_key = None

        # Motifs du catalogue (orientés) renommés dans le cadre de lecture du YAML :
        # {motif catalogue: motif YAML}, ex. {"GCA": "CAG"} pour ATXN1 (catalogue TGC, RC)
        self.motif_frame = {}

        # Samples pour ce TRID
        self.samples = {}

    def __repr__(self):
        return f"TRID({self.trid}, {len(self.samples)} samples)"
