
The config file lists paths to various files that are used by the workflow. This workflow only runs on hg38 data.

SCRIPT_NAME: The name of the compute job.
RESULTS_DIR: The output directory to store results.
INPUT_VCF: The input vcf file
SPLICEAI_PRECOMPUTED_ALL_VCF: The set of precomputed SpliceAI scores, downloaded from Illumina BaseSpace.
SPLICEAI_PRECOMPUTED_ALL_BASEPATH: The directory and base path for per chr precomputed SpliceAI scores.
REFERENCE_FASTA: The reference fasta file.
REFERENCE_FASTA_HG38: The hg38 reference fasta file. Likely the same as above.
REFERENCE_FASTA_HG19: The hg19 reference fasta file. This is used to liftover variants.
CHAIN_HG38_HG19: The chain file used to liftover variants from hg38 to hg19.
CLINVAR_VCF: ClinVar annotation as a vcf file (hg38).
HGMD_VCF: HGMD annotation as a vcf file (hg38).
GENES_BED_BASEPATH: The directory and base path for per chr hg38 bed files to be queried. Derived from SpliceAI annotation on their GitHub page.
BRANCHPOINT_MERCER_BED_HG38: The hg38 bed file annotating branchpoint regions in Mercer et al., 2015.
BRANCHPOINT_LABRANCHOR_BED_HG38: The hg38 bed file annotating branchpoint regions in Paggi et al., 2017 (LaBranchoR).
LCR_BED_HG38: The hg38 bed file annotating low complexity regions (LCR) used by gnomAD.
REPEATMASKER_BED_HG38: The hg38 bed file annotating RepeatMasker regions.
PAR_BED_HG38: The hg38 bed file annotating pseudoautosomal regions. Derived from Wikipedia.
GNOMAD2E_DIR: The directory containing gnomAD v2 exome vcfs, lifted over to hg38.
GNOMAD2G_DIR: The directory containing gnomAD v2 genome vcfs, lifted over to hg38.
GNOMAD3_DIR: The directory containing gnomAD v3 genome vcfs.
GNOMAD_CONSTRAINT: The path to the gnomAD v2 gene constraint file.
GNOMAD_COVERAGE: The path to the gnomAD v3 coverage file.
MANE_GRCH38: The path to the MANE canonical transcripts file.
GENES_CHD_ANNOTATION: The path to the file containing additional CHD gene annotations.
POOL_NAME: If the input vcf contains no samples (i.e. a pooled set of variants), a name can be provided here.


First execute the SpliceAI script with -c config.sh:
bash submit_spliceai.sh -c config.sh

Then summarize the annotated variants into a single R object:
bash submit_summarize_annotations.sh -c config.sh

If you want to only re-annotate filtered variants with VEP:
bash submit_reAnnotateVep.sh -c config.sh