SCRIPTS
bio_tools.py
Can be used to find gene ORF, transcript, and genomic (pre-splicing) sequences, and to align them (given CDS locations) in preparation for compute_features.precompute.
Also contains code for querying dbSNP, converting coordinates, etc.

compute_features.py
Uses a variety of tools to compute gene and variant features.
Will output wild type gene features (position-dependent conservation scores, domain. secondary structures, accessible surface area, etc) as a spreadsheet. 

parse_lovd2.py
Parses LOVD for combinations of disease-associated variants from patients, given a disease/gene. Then, queries dbSNP for benign (or at least non-pathogenic) variants to generate comparable control combinations. 

score_variants.py
Input is spreadsheet containing the variant combinations as rows and variants as columns, where 1 indicates variant contained in combination.
Utilizes gene spreadsheet, VEP, and other tools to compute variant scores, then combines multiple variants based on model.

combine_documents.py
Used to simplify analyses by reducing dimensions.

tools.py
Contains a variety of tools used in other scripts, including parsing variant specifics from a string, finding the associated AA variant for an NT variant, aligning multiple sequences, etc.

DEPENDENCIES 
Clustal Omega
NetSurfP2
netphos
netNglyc
Nupack
mFold
KineFold
remuRNA
PLMC
Python libraries: pandas, Biopython, scikit-learn, numpy, scipy, numba, networkx, prody, BeautifulSoup, intervaltree

RARE CODON ENRICHMENT CALCULATION (rc_enrichment)
Originally from Jacobs and Shakhnovich, 2017 (https://faculty.chemistry.harvard.edu/shakhnovich/software/coarse-grained-co-translational-folding-analysis), but optimized here to improve Poisson binomial calculation. 
Requires PLMC installation.

EVMUTATION (EVmutation)
Originally from Marks lab (Hopf et. al., 2017) (https://github.com/debbiemarkslab/EVmutation). Given here with a wrapper to realign with the focus sequence.

HOW TO USE
We strongly recommend having a disk mounted in RAM, as many intermediate programs use disk writes, which can reduce the lifespan of SSDs. This variable is specified in each .py file as "ram_disk".
Note that compute_features.precompute only needs to be run once for each gene. For genes involved in multiple diseases (i.e. VWF), this can speed up runtime.
Automated querying of Conseq does not work, and a GENE_consurf.tsv file must be provided in the RAM disk to find Conseq scores.

for (GENE, DISEASE_ID, REFSEQ_ACCESSION) in GENESET:
	bio_tools.precompute(GENE, outdir=...) #to write out the aligned ORF, transcript, and genomic sequences for the gene to the outdirectory as GENE.fasta
	seqs = compute_features.parse_seqs(outdir, GENE) #to load these sequences
	genedata, ids = compute_features.precompute(gene, seqs, path=/path/to/RAM/disk) #to precompute the gene conservation, domain, secondary structure, accessible surface area, etc
	parse_lovd2.pipeline(GENE, DISEASE_ID, REFSEQ_ACCESSION) #reads all patient variant combinations from LOVD for certain disease (assuming correct RefSeq accession)
score_variants.pipeline([(GENE1, OUTPUT_DIR1, LOVD_OUTPUT1), ..., (GENEn, OUTPUT_DIRn, LOVD_OUTPUTn)]) #computes gene features for all variant combinations in all genes/diseases
