# TGV — TRGT Global Viewer

![Pipeline TGV](assets/workflow_simple.png)

<details open>
  <summary><b>🇫🇷 Version Française (Cliquez pour replier)</b></summary>
  <br>

**TGV (TRGT Global Viewer)** est un outil de visualisation clinique conçu pour le CHU de Nîmes. Il simplifie l'analyse, le contrôle qualité et l'interprétation des répétitions en tandem issues du workflow **TRGT (PacBio SMRT Link)**.

---

### 🧬 Contexte Clinique & Vision "Zéro Désarchivage"

* **Le défi diagnostique** : L'analyse des expansions de répétitions en tandem (ex: ataxies) repose sur le séquençage HiFi de haute précision à longues lectures (séquenceur **PacBio Vega**). Si l'outil TRGT offre un profilage génomique puissant, les données brutes générées sont denses, éparpillées et complexes à manipuler.
* **La vision "Zéro désarchivage manuel" (Ergonomie)** : Au quotidien, manipuler et décompresser manuellement des dizaines d'archives ZIP volumineuses (fichiers BAM de plusieurs centaines de mégaoctets, graphiques d'allèles ou de méthylation) est une tâche lourde. **TGV résout ce problème en lisant et en traitant toutes les données directement à la volée en arrière-plan, sans aucune décompression manuelle préalable sur le disque dur de l'utilisateur.**

---

### 📋 Caractéristiques principales

* **Interface graphique (GUI) intuitive** :  Chargement des données, filtrage par patient ou locus (TRID), et surcharge manuelle tracée des classifications et génotypes (les seuils cliniques se modifient dans `clinical_thresholds.yaml`).
* **Cochage automatique par panel** : Intègre un système de boutons configurables permettant de sélectionner automatiquement en un clic des listes de gènes d'intérêt (panels in silico) définies par l'utilisateur.
* **Visualisation d'alignements (igv.js)** :  Extraction automatique des bams (spanning_BAM / repeat_reads) et ouverture d'une session IGV locale dans votre navigateur pour une inspection immédiate.
* **Affichage de graphiques TRGT (SVG)** : Rendu direct des profils d'allèles et de méthylation générés par l'outil TRVZ (TRGT), sans extraction manuelle.
* **Rapport de contrôle qualité global (QC)** : Rapport HTML généré à la volée pour visualiser la qualité d'enrichissement du run. Note : QC produit via le module tgv_inputs_builder.py.
* **Traçabilité par fichier de logs** : Chaque analyse génère un journal structuré et horodaté dans logs/ (versions, entrées, étapes du pipeline et alertes).
* **Fichiers temporaires confinés** : Les fichiers extraits (BAM pour IGV, graphiques, exports HTML, rapport QC) sont écrits dans un répertoire temporaire propre à la session (`tgv_*`), supprimé à la fermeture de l'application. Ils ne sont servis au navigateur que par un serveur local restreint (127.0.0.1, jeton aléatoire, liste blanche de fichiers).
* **Léger et portable** : Développé sans dépendances lourdes (pas de Pandas, NumPy ou Jinja2). Disponible en exécutable autonome (Windows) ou script léger (Linux/macOS).

---

### 🚀 Comment l'utiliser (Usage)

L'outil **TGV** s'adapte à votre environnement via deux modes d'exécution :

#### Option A : Sous Windows (Exécutable autonome)
Destiné aux cliniciens et biologistes sur poste de travail Windows.
1. Récupérez l'exécutable autonome **`TGV.exe`** produit par le workflow GitHub Actions *Build Windows EXE* (onglet *Actions*, artefact `TGV-exe`). Aucune release n'est publiée à ce jour.
2. Double-cliquez sur l'exécutable pour lancer l'application.
*Aucune installation de Python ou de bibliothèque n'est requise.*

#### Option B : Sous Linux / macOS (Ligne de commande)
Destiné aux bio-informaticiens ou pour une utilisation sur serveur de calcul.

1. Installez les dépendances (versions épinglées) :
   ```bash
   pip install -r requirements.txt            # interface graphique et CLI
   pip install -r requirements-builder.txt    # tgv_inputs_builder.py (json5, matplotlib)
   ```
2. Lancez l'application :
   ```bash
   python main.py
   ```

#### Option C : Ligne de commande (CLI, intégration cron / Slurm)
```bash
python tgv_cli.py --zip RUN-trgt_vcfs.zip --out export.tsv            # panel « Ataxie »
python tgv_cli.py --zip RUN-trgt_vcfs.zip --out export.tsv --all-loci # tous les loci
```
* La CLI ne nécessite ni tkinter ni affichage graphique.
* **Codes retour** : `0` succès ; `1` au moins un échantillon en échec (les autres sont écrits) ou aucune ligne produite ; `2` configuration ou entrée invalide (rien n'est écrit).
* **Relance idempotente** : les lignes existantes d'un même couple (`run_id`, `sample_id`) sont remplacées ; `run_id` vaut par défaut le préfixe du ZIP.
* **Provenance** : chaque ligne porte la version de TRGT et le catalogue (lus dans l'en-tête du VCF), la version et le commit de TGV, l'empreinte SHA-256 de `clinical_thresholds.yaml` et la référence des seuils (`source`).

#### Reproduire l'analyse avec le catalogue public de TRGT
```bash
trgt genotype --genome GRCh38.fa --repeats pathogenic_repeats.hg38.bed \
    --reads S1.bam --output-prefix S1.trgt --preset targeted --karyotype XY
bcftools sort -Oz -o S1.trgt.sorted.vcf.gz S1.trgt.vcf.gz
gunzip -c S1.trgt.sorted.vcf.gz > S1.trgt.sorted.vcf
zip RUN-trgt_vcfs.zip S1.trgt.sorted.vcf
python tgv_cli.py --zip RUN-trgt_vcfs.zip --out export.tsv --all-loci
```
Le catalogue public ([`pathogenic_repeats.hg38.bed`](https://github.com/PacificBiosciences/trgt/blob/main/repeats/pathogenic_repeats.hg38.bed)) utilise des identifiants de gènes (`ATXN1`, `TBP`…) : `configs/trid_aliases.yaml` les relie aux blocs de `clinical_thresholds.yaml`. Des jeux de données PacBio publics sont disponibles, par exemple sur [pacb.com/vega-targeted-datasets](https://www.pacb.com/vega-targeted-datasets/).

---

### 📂 Spécifications des Entrées (Inputs)

Pour fonctionner de manière transparente et sans désarchivage préalable, **TGV** s'appuie sur une structure d'archives standardisée **au niveau du Run** (les fichiers d'échantillons individuels sont stockés à l'intérieur d'archives globales du run) :

#### 1. Fichiers généraux de l'analyse (Sélection manuelle sur l'IHM)
* **Archive de VCFs TRGT (`{ID_RUN}-trgt_vcfs.zip`)** : L'archive ZIP contenant l'ensemble des fichiers `.trgt.vcf` du run (un fichier VCF par échantillon).
* **Génome de référence (Optionnel)** : Le fichier de référence génomique au format `.fa` ou `.fasta` accompagné de son fichier d'index `.fai` (ex : `hg38.fa` et `hg38.fa.fai`).

#### 2. Archives globales de Run détectées automatiquement (Sister ZIPs)
TGV détecte automatiquement les archives associées présentes dans le même répertoire, sous réserve qu'elles partagent exactement le même préfixe de run (`{ID_RUN}-`) :
* **Alignements (BAM)** :
  * `{ID_RUN}-spanning_BAM.zip` : Contient les fichiers BAM/BAI de type *spanning* pour tous les patients du run.
  * `{ID_RUN}-repeat_reads.zip` : Contient les fichiers BAM/BAI de type *mapped* pour tous les patients du run.
* **Graphiques TRGT (SVG)** :
  * `{ID_RUN}-trgt_motifs_allele.zip` : Archives des profils de tailles des motifs d'allèles de tous les patients.
  * `{ID_RUN}-trgt_motifs_waterfall.zip` : Archives des profils de reads *waterfall* de tous les patients.
  * `{ID_RUN}-trgt_meth_allele.zip` : Archives de méthylation allèle-spécifique de tous les patients.
  * `{ID_RUN}-trgt_meth_waterfall.zip` : Archives de profils de reads *waterfall* de méthylation de tous les patients.

*Fonctionnement interne : Lors de la sélection d'un patient et d'un locus, TGV ouvre l'archive globale du run correspondante, y recherche le fichier du patient **par correspondance exacte** sur les conventions de `tgv_inputs_builder.py` (`{patient}.trgt.spanning.sorted.bam`, `{patient}.repeat_reads.bam`, `{patient}_{catégorie}.trvz_alleles.zip`), l'extrait dans le répertoire temporaire de la session, puis le supprime à la fermeture. Si plusieurs fichiers correspondent, aucun n'est ouvert. Le nom réel du fichier est affiché dans IGV.*

#### 3. Rapport de Contrôle Qualité (Module exclusif TGV)
Ce rapport, qui n'est pas généré nativement par SMRT Link/TRGT, est **produit spécifiquement par votre script tgv_inputs_builder.py** à partir des données brutes du run. Il permet une supervision globale que le workflow standard ne propose pas.
L'archive **{id}-QC.zip** (générée par le builder) contient :
* **Rapport d'enrichissement** : target_enrichment_puretarget.report.json
* **Synthèse échantillons** : sample_summary.csv
* **Couverture par cible** : target_cov_by_sample.csv
* **Visualisations** : Graphiques PNG (ex: sample_coverage_boxplot-0.png, read_categories.png).

---

### ⚙️ Configuration & Personnalisation

TGV est configurable grâce aux fichiers du répertoire **`configs/`** :

* **`clinical_thresholds.yaml`** : Fichier de référence clinique. Il définit, pour chaque maladie/locus (TRID), les plages de répétitions permettant de classer les allèles ainsi que l'orientation du brin (directe ou reverse-complement). Il est **validé au chargement** : une erreur (label absent de `label_priority`, borne invalide, orientation manquante…) bloque l'analyse ; les trous et chevauchements de plages sont signalés.
* **`buttons_panel.yaml`** : Panels de loci sélectionnables en un clic (ex : « Ataxie »), également utilisés par la CLI.
* **`trid_aliases.yaml`** : Correspondance entre les identifiants d'un autre catalogue TRGT (ex : catalogue public) et les blocs de `clinical_thresholds.yaml`.
* **`trgt_params.json5`** : Paramètres de `trgt genotype` et `trgt plot` utilisés par `tgv_inputs_builder.py` (obligatoire : sans lui, le builder s'arrête plutôt que de lancer TRGT avec le preset par défaut `wgs`).

#### Statuts particuliers
* **`unclassified`** : allèle appelé mais non classable (valeur hors de toute plage du YAML, plages qui se chevauchent, ou motifs TRGT absents des groupes cliniques). TGV ne le classe jamais « normal » par défaut ; le locus reste dans l'export avec l'explication en commentaire.
* **`no_call`** : allèle non appelé par TRGT ; **`-`** : second allèle absent d'un locus haploïde (ex. chrX, caryotype XY).
* Les deux allèles sont affichés dans l'ordre croissant du génotype clinique (allèle 1 = plus petit).

---

### 🛠️ Organisation des fichiers de logs

Pour assurer une traçabilité complète de vos analyses, **TGV** génère automatiquement des journaux horodatés dans le sous-dossier `logs/` :

*   **Pour l'interface graphique** :
    `tgv.ANNEEMOISJOURHEUREMINUTESECONDE.log` (ex: `tgv.20260612145002.log`)
*   **Pour la préparation des données (Builder)** :
    `tgv_input_builder.ANNEEMOISJOURHEUREMINUTESECONDE.log` (ex: `tgv_input_builder.20260922135222.log`)
*   **Pour la CLI** : journal sur la sortie console (à rediriger par l'ordonnanceur).

Ces fichiers contiennent la version et le commit de TGV, l'empreinte SHA-256 de `clinical_thresholds.yaml`, les avertissements de configuration et le journal d'audit des surcharges manuelles (utilisateur, patient, run, valeur automatique d'origine).

---

### ⚠️ Limites connues

* Les règles cliniques suivantes restent à arbitrer par les biologistes et ne sont pas modifiées par le code : SCA1 (toute interruption rend un allèle de 39–44 unités « normal », GeneReviews ne retient que les interruptions CAT), FGF14 (expansions GAA non pures), FXN (`protective_motifs` non utilisé), écarts de seuils avec GeneReviews.
* La règle de profondeur (seuil par allèle) pénalise les homozygotes, dont les lectures sont réparties entre deux allèles identiques.
* `clinical_thresholds.yaml` comporte des trous (RFC1 < 200, SCA7 20–27, SCA36 15–649…) et un chevauchement (FMR1 = 200) : les valeurs concernées sont classées `unclassified`.
* PySimpleGUI-4-foss est un miroir figé de PySimpleGUI 4, sans mises à jour.

---

### 🧪 Tests
```bash
pip install -r requirements-dev.txt
ruff check .
pytest
```
Les tests (VCF synthétiques au format TRGT 5.x) couvrent le parsing, la classification, la CLI, le serveur local et le builder ; ils tournent sans interface graphique et sont exécutés par la CI (`.github/workflows/tests.yaml`).

</details>

<br>

<details>
  <summary><b>🇬🇧 English Version</b></summary>
  <br>

**TGV (TRGT Global Viewer)** is a clinical visualization tool designed for Nîmes University Hospital. It streamlines analysis, quality control, and interpretation of tandem repeats from the **TRGT (PacBio SMRT Link)** workflow.

---

### 🧬 Clinical Context & "Zero Extraction" Philosophy

* **The diagnostic challenge**: Interpreting expansions of tandem repeats relies on high-precision HiFi sequencing (**PacBio Vega** sequencer). While TRGT provides powerful genomic profiling, raw outputs are dense and complex to handle.
* **The "Zero manual extraction" philosophy (Ergonomics)**: Handling dozens of heavy ZIP archives daily (BAM files, allele plots, methylation data) is tedious. **TGV processes all data on-the-fly in the background, eliminating the need for manual unarchiving on your local drive.**

---

### 📋 Main Features

* **Intuitive GUI**: Easily load data, filter patients or loci (TRIDs), and apply traced manual overrides of classifications and genotypes (clinical thresholds are edited in `clinical_thresholds.yaml`).
* **Automation Modules**: 
    * **`tgv_inputs_builder.py`**: Automates archive structuring and QC report generation.
    * **`tgv_cli.py`**: Command-line interface for batch processing and pipeline integration.
* **Alignment Visualization (igv.js)**: Automated BAM extraction, served by a restricted local HTTP server (127.0.0.1, random token, file whitelist, Range requests).
* **TRGT Graphics Display (SVG)**: Direct rendering of allele/methylation plots (TRVZ tool) with no manual extraction.
* **Global Run QC**: On-the-fly HTML quality report (generated via `tgv_inputs_builder.py`).
* **Traceability & Logs**: Structured, timestamped logs for both the GUI and the Builder in `logs/`.
* **Confined temporary files**: One temporary directory per session (`tgv_*`), removed on exit.
* **Lightweight & Portable**: No heavy dependencies (no Pandas/NumPy). Portable executable (Windows) or simple script (Linux/macOS).

---

### 🚀 How to use (Usage)

**TGV** adapts to your work environment through two execution modes:

#### Option A: Windows (Standalone executable)
Aimed at clinicians and biologists on Windows workstations.
1. Get the standalone **`TGV.exe`** from the *Build Windows EXE* GitHub Actions workflow (artifact `TGV-exe`); no release has been published yet.
2. Double-click to launch.

#### Option B: Linux / macOS (Command-line usage)
Aimed at bioinformaticians or server environment usage.
1. Install pinned dependencies: `pip install -r requirements.txt` (and `requirements-builder.txt` for the builder)
2. GUI Mode: `python main.py` | CLI Mode: `python tgv_cli.py --help`
3. CLI exit codes: `0` success, `1` sample failure or no rows, `2` invalid configuration/input. Re-running replaces the rows of the same (`run_id`, `sample_id`). Each row carries TRGT version and catalog (VCF header), TGV version/commit and the SHA-256 of `clinical_thresholds.yaml`.

---

### 📂 Input Specifications

**TGV** relies on standardized Run-level ZIP archives:
* **Manual Inputs (GUI)**: TRGT VCF archive (`{ID_RUN}-trgt_vcfs.zip`) and optional Reference Genome.
* **Auto-detected Sister ZIPs**: Alignment files (`spanning_BAM`, `repeat_reads`) and TRGT plot archives (`motifs_allele`, `motifs_waterfall`, `meth_allele`, `meth_waterfall`).
* **QC Report (Exclusive TGV Module)**: Archive `{id}-QC.zip` (produced by `tgv_inputs_builder.py`) containing enrichment reports, summaries, and coverage plots.

---

### ⚙️ Configuration & Customization

* **`clinical_thresholds.yaml`**: Defines clinical repeat ranges and motif strand orientation; validated on load (errors block the analysis, gaps/overlaps are reported and yield `unclassified`).
* **`buttons_panel.yaml`**: Allows customizing the GUI with dynamic buttons for specific gene panels (e.g., "Ataxias").
* **`trid_aliases.yaml`**: Maps TRIDs of other TRGT catalogs (e.g., the public catalog) to `clinical_thresholds.yaml` blocks.
* **`trgt_params.json5`**: TRGT parameters used by `tgv_inputs_builder.py` (required).

Allele statuses: `unclassified` (called but not classifiable — never defaulted to "normal"), `no_call` (not called by TRGT), `-` (absent second allele at a haploid locus). Alleles are ordered by increasing clinical genotype.

---

### 🛠️ Logs Directory Organization

For full auditability, **TGV** generates timestamped logs in the `logs/` subdirectory:
*   **GUI**: `tgv.YYYYMMDDHHMMSS.log`
*   **Data Preparation (Builder)**: `tgv_input_builder.YYYYMMDDHHMMSS.log`
*   **CLI**: console output.

Known limitations (clinical rules for SCA1/FGF14/FXN, depth rule for homozygotes, gaps in `clinical_thresholds.yaml`, frozen PySimpleGUI-4-foss) are listed in the French section. Tests: `pip install -r requirements-dev.txt && ruff check . && pytest`.

</details>

<br>

<details>
  <summary><b>📚 Références & Outils tiers / References & Third-Party Tools</b></summary>
  <br>

#### 🇫🇷 Références & Outils tiers
Si vous utilisez **TGV** dans le cadre de vos travaux cliniques ou de recherche, veuillez citer les outils sous-jacents :

* **TRGT / TRVZ** :
  > Dolzhenko, E., English, A., Dashnow, H., *et al.* Characterization and visualization of tandem repeats at genome scale. *Nat Biotechnol* (2024). [https://doi.org/10.1038/s41587-023-02057-3](https://doi.org/10.1038/s41587-023-02057-3)
* **igv.js** :
  > Robinson, J. T., Thorvaldsdóttir, H., Turner, D., & Mesirov, J. P. igv.js: an embeddable JavaScript implementation of the Integrative Genomics Viewer (IGV). *Bioinformatics*, 39(1), btac830 (2023). [https://doi.org/10.1093/bioinformatics/btac830](https://doi.org/10.1093/bioinformatics/btac830)

---

#### 🇬🇧 References & Third-Party Tools
If you use **TGV** for clinical work or scientific publications, please acknowledge the underlying tools:

* **TRGT / TRVZ**:
  > Dolzhenko, E., English, A., Dashnow, H., *et al.* Characterization and visualization of tandem repeats at genome scale. *Nat Biotechnol* (2024). [https://doi.org/10.1038/s41587-023-02057-3](https://doi.org/10.1038/s41587-023-02057-3)
* **igv.js**:
  > Robinson, J. T., Thorvaldsdóttir, H., Turner, D., & Mesirov, J. P. igv.js: an embeddable JavaScript implementation of the Integrative Genomics Viewer (IGV). *Bioinformatics*, 39(1), btac830 (2023). [https://doi.org/10.1093/bioinformatics/btac830](https://doi.org/10.1093/bioinformatics/btac830)

</details>

<br>

<details>
  <summary><b>📝 Licence & Auteurs / License & Authors</b></summary>
  <br>

* **Auteur principal / Main Author** : Corentin Marco (CHU de Nîmes)
* **Licence / License** : [MIT](LICENSE).

</details>
