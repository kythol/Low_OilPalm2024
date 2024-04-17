# Low_OilPalm2024 (G3)

This repository contains supplementary materials (script) used for model classification in _Low et al 2024 (G3) Chromosome-scale _Elaeis guineensis_ and _E. oleifera_ assemblies: Comparative genomics of oil palm and other Arecaceae_

Usage:
Rscript to_publish_P11_overall_models_script.r input/testset_seqping.gff3 \ # MAKER predictions
    input/testset_mikado.gff3 \ #mikado predictions
    input/testset_TSS.osc \ # TSS information
    input/testset_mikado_proper_ATG.txt \ # gene models list with good ATG
    input/testset_gtf \ # folder with expression GTFs
    input/testset_maker_Liliopsida.txt \ #BLAST results for maker
    input/testset_mikado_Liliopsida.txt \ # BLAST results for mikado
    input/testset_overlaps.txt # list of overlapping models

Organisation of output: Two folders are created: tmp and results. The main output (after classification and representative steps) are saved into results/, the other 6 intermediate files (two CAGE files, file after integration, after expression data added, after blast data added, after overlaps added) are saved in tmp/ folder.

Summary:
Classification: First, we analyse maker and mikado predictions separately: 1. ATG coordinates are identified for each isoform based on the first CDS location (this is NOT a filtering step, just adding information); 2. Association of CAGE tags with isoforms: the inclusion criterion is that the CAGE tag is located within -3000 to -10 from the ATG for the gene isoform. Only one tag should be used per gene, and if the tag is associated for the gene, it can't be used for other genes, but can be used for multiple isoforms of one gene. In case there are several tags passing the criteria, the closest tag to ATG is picked. This is the part where we additionally check the mikado isoforms for a start/stop codon (this is NOT a filtering step, just adding information); 3. Model integration: the goal is to integrate mikado and maker results into one set of predictions. First, in case of obviously identical models that are shared by two datasets, the mikado model is removed. This is decided based on the CDS coordinates of the models. If the mRNA coordinates are identical but the CDS coordinates are different, we include both and assign them to the same Parent ID (maker's Parent ID). (if they are identical, the CAGE tag should be associated to the maker model, and not the mikado model. If they are similar but not identical, e.g. the ATGs are different, there can be a different CAGE tag for mikado model. Overlapping models are not considered yet!) 4. Expression data: for each GTF file, the TPM counts for transcripts were extracted, and the transcript coordinates are compared with the gene models. The transcript should overlap the gene models with at least 80% coverage of it's length; these ratios are reported in the final file along with TPM counts. The data is then added to the integrated gene model list. 5. BLAST: blast search performed with the following parameters: -evalue 0.00001 -maxtargetseqs 5 -word_size 4 -threshold 21. The data is added to the integrated list if it is present for the isoform. 6. Classification of models: for classification, the evidence is converted into binary format: if present it is a "+", if not present it's a ".". Expression data is considered to be present if the overlap of the transcript with the gene models is >= 80% and the TPM is >= 1. The CAGE data is considered to be present if there is a tag associated with the gene model. BLAST data is considered in the same way. Based on these criteria, the following classes are identified: Class 1: models with CAGE, BLAST & expression (3 out of 3 evidences) Class 2: models with CAGE & BLAST or CGE and expression (2 out of 3) Class 3: models with BLAST & expression (2 out of 3) Class 4: models with only BLAST (1 out of 3) Class 5: models with either CAGE or expression (1 out of 3) Class 6: no evidence

Representatives: Takes the output of the classification script. First, we correct the Parent for the models that were predicted using different tools, but are actually part of the same gene. We identify the representatives only for multi-isoform genes; the single-isoform genes are automatically labelled as representatives. The main criteria to consider the representative gene is the category given during classification (priority to the categories 1, 2, 3 etc). In case of several genes meeting the first criteria, the models are tested for the BLAST e-value (for categories 1, 2 and 4), followed by the expected CAGE tag density (based on sorgum). Then, for the number of tissues expressing this isoform, followed by max TPM counts, max CDS length. In case we still have several candidates, the first one assigned by the prediction tool is used (_RA for maker and .1 for mikado). If maker and mikado predictions are present, we prioritise maker prediction.
