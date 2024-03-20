#!/usr/bin/env Rscript

#################################################################
# Developer: Robert Lesurf
# Arguments: config file containing various directories and paths
# Take annotated SpliceAI vcf files and further summarize
#################################################################


##### Read in the config file #####
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=1)
  stop("A config file must be suppled as a command line argument")
CONFIG = read.table(args[1], sep="=")
for(i in 1:nrow(CONFIG))
  assign(as.character(CONFIG[i,1]), as.character(CONFIG[i,2]))

# set defaults as needed
if(!exists("RESULTS_DIR"))
  stop("The RESULTS_DIR variable is missing in the config file.")



##### Check that all required modules have been loaded #####

err <- tryCatch(
    {
    bcftoolsVersion <- system("bcftools --version", intern=TRUE)
    }, warning = function(w) {
        return("WARNING")
    }, error = function(e) {
        return("ERROR")
    }
)
if(is.null(err) || err[1] == "ERROR")
    stop("You must load the following modules before running this script: bcftools/1.9")

err <- tryCatch(
    {
    javaVersion <- system("gatk --version", intern=TRUE)
    }, warning = function(w) {
        return("WARNING")
    }, error = function(e) {
        return("ERROR")
    }
)
if(is.null(err) || err[1] == "ERROR")
    stop("You must load the following modules before running this script: gatk/4.2.2.0")



##### Helper functions #####

get_vcf_variants = function(vcf) {
  variants = system(paste0("bcftools query -f '%CHROM %POS %REF %ALT\n' ", vcf), intern=TRUE)
  #variant_keys = gsub(" ", "_", variants)
  variant_kcpra = t(sapply(variants, function(x){
    x_parsed = unlist(strsplit(x, "[ ]"))
    return(c("Key"=gsub(" ","_",x), "CHROM"=x_parsed[1], "POS"=x_parsed[2], "REF"=x_parsed[3], "ALT"=x_parsed[4]))
  }))
  rownames(variant_kcpra) = NULL
  return(as.data.frame(variant_kcpra))
}


get_vcf_spliceai = function(vcf) {
  variant_keys = get_vcf_variants(vcf)$Key
  if(length(variant_keys) == 0)
    return()

  # A single variant can have multiple genes/interpretations, separated by a comma (no spaces)
  # ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
  spliceai_format = c("SpliceAI_ALLELE", "SpliceAI_SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL")

  info_spliceai = system(paste0("bcftools query -f '%INFO/SpliceAI\n' ", vcf), intern=TRUE)
  info_spliceai_split = lapply(info_spliceai, function(x){unlist(strsplit(x,"[,]"))})
  spliceai_df = matrix(nrow=sum(sapply(info_spliceai_split,length)), ncol=1+10, dimnames=list(NULL,c("Key",spliceai_format)))
  pointer = 1
  for(i in 1:length(info_spliceai_split)) {
    for(j in 1:length(info_spliceai_split[[i]])) {
      ij_parsed = unlist(strsplit(info_spliceai_split[[i]][j], "[|]"))
      ij_parsed[which(ij_parsed == "")] = NA
      spliceai_df[pointer, "Key"] = variant_keys[i]
      spliceai_df[pointer, 2:(1+length(ij_parsed))] = ij_parsed
      pointer = pointer + 1
    }
  }
  spliceai_df = as.data.frame(unique(spliceai_df))
  return(spliceai_df)
}


get_vcf_vepCsq = function(vcf, mane=MANE_GRCH38) {
  variant_keys = get_vcf_variants(vcf)$Key

  # Parse the CSQ header field
  header = system(paste0("bcftools view -h ", vcf), intern=TRUE)
  csq_ind = which(sapply(header, substr, 0, 15) == "##INFO=<ID=CSQ,")
  if(length(csq_ind) != 1)
    stop("Cannot find a single CSQ info filed in the header")
  csq_format = unlist(strsplit(header[csq_ind], "[:]"))[2]
  csq_format = substr(csq_format, 2, nchar(csq_format) - 2) # starts with a blank space and ends with ">
  csq_format = unlist(strsplit(csq_format, "[|]"))

  # Extract and label the variant CSQ data. A single variant can have multiple comma-separated CSQ annotations
  info_csq = system(paste0("bcftools query -f '%INFO/CSQ\n' ", vcf), intern=TRUE)
  info_csq_split = lapply(info_csq, function(x){unlist(strsplit(x,"[,]"))})
  csq_df = matrix(nrow=sum(sapply(info_csq_split,length)), ncol=1+length(csq_format), dimnames=list(NULL,c("Key",csq_format)))
  pointer = 1
  for(i in 1:length(info_csq_split)) {
    for(j in 1:length(info_csq_split[[i]])) {
      ij_parsed = unlist(strsplit(info_csq_split[[i]][j], "[|]"))
      ij_parsed[which(ij_parsed == "")] = NA
      csq_df[pointer, "Key"] = variant_keys[i]
      csq_df[pointer, 2:(1+length(ij_parsed))] = ij_parsed
      pointer = pointer + 1
    }
  }
  csq_df = as.data.frame(unique(csq_df))

  # Annotate MANE transcripts
  maneSummary = read.delim(MANE_GRCH38)
  csqFeature_parsed = sapply(csq_df$Feature, function(x){unlist(strsplit(x,"[.]"))[1]})
  maneRefseq_parsed = sapply(as.character(maneSummary$RefSeq_nuc), function(x){unlist(strsplit(x,"[.]"))[1]})
  maneEnst_parsed   = sapply(as.character(maneSummary$Ensembl_nuc), function(x){unlist(strsplit(x,"[.]"))[1]})
  maneInd = sapply(csqFeature_parsed, function(x){
    xInd = union(which(maneRefseq_parsed==x), which(maneEnst_parsed==x))
    if(length(xInd)==1)
      return(xInd)
    return(NA)
  })
  colnames(maneSummary)[which(colnames(maneSummary)=="symbol")] = "MANE_Symbol"
  colnames(maneSummary)[which(colnames(maneSummary)=="X.NCBI_GeneID")] = "MANE_EntrezID"
  colnames(maneSummary)[which(colnames(maneSummary)=="Ensembl_Gene")] = "MANE_ENSG"
  colnames(maneSummary)[which(colnames(maneSummary)=="HGNC_ID")] = "MANE_HGNC"
  colnames(maneSummary)[which(colnames(maneSummary)=="RefSeq_nuc")] = "MANE_RefSeq"
  colnames(maneSummary)[which(colnames(maneSummary)=="Ensembl_nuc")] = "MANE_ENST"

  # For each key: if there are mane and non-mane transcripts, keep only the mane transcripts
  csqmane = cbind(csq_df, maneSummary[maneInd,c("MANE_status","MANE_Symbol","MANE_EntrezID","MANE_ENSG","MANE_HGNC","MANE_RefSeq","MANE_ENST")])
  csqmane_ind = unique(unlist(lapply(unique(variant_keys), function(x){
    x_ind = which(csqmane == x)
    x_ind_mane = x_ind[which(!is.na(csqmane$MANE_status[x_ind]))]
    if(length(x_ind_mane) > 0)
      return(x_ind_mane)
    return(x_ind)
  })))
  csqmane_reduced = csqmane[csqmane_ind, ]
  csqmane_reduced = droplevels(csqmane_reduced)

  return(csqmane_reduced)
}


get_variant_liftover = function(vcf, scratch_dir, chain=CHAIN_HG38_HG19, fasta=REFERENCE_FASTA_HG19) {
  if(!dir.exists(scratch_dir))
    dir.create(scratch_dir)
  currenttime = format(Sys.time(),"%Y-%m-%d-%H-%M-%S")
  vcf_out = paste0(scratch_dir, "/liftover_", currenttime, "_", basename(vcf))
  vcf_reject = paste0(scratch_dir, "/liftoverReject_", currenttime, "_", basename(vcf))

  cmd_liftover = "java -Xmx5g -jar /hpf/tools/centos7/gatk/4.2.2.0/gatk-package-4.2.2.0-local.jar LiftoverVcf"
  cmd_liftover = paste(cmd_liftover, "-I",vcf, "-O",vcf_out, "--REJECT",vcf_reject)
  cmd_liftover = paste(cmd_liftover, "--CHAIN",chain, "-R",fasta)
  cmd_liftover = paste(cmd_liftover, "--WRITE_ORIGINAL_POSITION --MAX_RECORDS_IN_RAM 50000")
  system(cmd_liftover)

  variant_in_keys = get_vcf_variants(vcf)$Key
  positions_in = system(paste0("bcftools query -f '%CHROM %POS\n' ", vcf), intern=TRUE)
  refalt_in = system(paste0("bcftools query -f '%REF %ALT\n' ", vcf), intern=TRUE)
  variant_out_keys = get_vcf_variants(vcf_out)$Key
  positions_out_original = system(paste0("bcftools query -f '%INFO/OriginalContig %INFO/OriginalStart\n' ", vcf_out), intern=TRUE)
  refalt_out = system(paste0("bcftools query -f '%REF %ALT\n' ", vcf_out), intern=TRUE)

  file.remove(vcf_out)
  file.remove(vcf_reject)

  ind_out = rep(NA, length(variant_in_keys))
  for(i in 1:length(ind_out)) {
    ind_tmp = unique(which(positions_out_original == positions_in[i]))
    if(length(ind_tmp) > 1)
      ind_tmp = intersect(ind_tmp, which(refalt_out == refalt_in[i]))
    if(length(ind_tmp) == 1)
      ind_out[i] = ind_tmp
  }
    
  liftover_df = data.frame("Key"=variant_in_keys, "Key_LiftOver"=variant_out_keys[ind_out])
  rownames(liftover_df) = NULL
  liftover_df = unique(liftover_df)
  return(liftover_df)
}


get_variant_bed = function(vcf, bed, label="") {
  if(label == "")
    stop(paste("You must provide a label for the bed file:", bed))
  variant_keys = get_vcf_variants(vcf)$Key
  variant_keys_bed = system(paste0("bcftools query -f '%CHROM %POS %REF %ALT\n' -R ", bed, " ", vcf), intern=TRUE)
  variant_keys_bed = gsub(" ", "_", variant_keys_bed)
  variant_bed_df = data.frame("Key"=variant_keys, "Label"=sapply(variant_keys, function(x){x %in% variant_keys_bed}))
  colnames(variant_bed_df) = c("Key", label)
  rownames(variant_bed_df) = NULL
  return(variant_bed_df)
}


get_variant_altSamples = function(vcf) {
  vcf_samples = system(paste0("bcftools query -l ", vcf), intern=TRUE)
  if(length(vcf_samples) > 0) {
    variantKeysSamplesGt = system(paste0("bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE=%GT ]\n' -i 'GT=\"alt\"' ", vcf), intern=TRUE)
    variantKeysSamplesGt_parsed = t(sapply(variantKeysSamplesGt, function(x){unlist(strsplit(x,"\t"))}))
    rownames(variantKeysSamplesGt_parsed) = NULL
    variantKeysSamplesGt_df = data.frame("Key"=apply(variantKeysSamplesGt_parsed[,1:4], 1, paste, collapse="_"), "Samples"=variantKeysSamplesGt_parsed[,5])
    return(unique(variantKeysSamplesGt_df))
  } else {
    variantKeys = system(paste0("bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ", vcf), intern=TRUE)
    variantKeys_parsed = t(sapply(variantKeys, function(x){unlist(strsplit(x,"\t"))}))
    rownames(variantKeys_parsed) = NULL
    variantKeysPooled_df = data.frame("Key"=apply(variantKeys_parsed[,1:4], 1, paste, collapse="_"), "Samples"=paste0(POOL_NAME,"=1"))
    return(unique(variantKeysPooled_df))
  }
}


get_gene_gnomadConstraint = function(enst=NULL, symbol=NULL) {
  # Note that some gene symbols have 2 entries in the gnomAD constraint file. If there are multiple matches, nothing gets returned.

  query_length = max(length(enst), length(symbol))
  if(!all(c(length(enst), length(symbol)) %in% c(0, query_length)))
    stop("Input query variables have different lengths")

  # Merge input queries into a single data.frame, and remove trailing .version labels for ENST transcripts
  query_df = data.frame(
    "ENST_split" = sapply(1:query_length, function(x){ifelse(length(enst)==0, NA, unlist(strsplit(as.character(enst[x]),"[.]"))[1])}),
    "Symbol" = sapply(1:query_length, function(x){ifelse(length(symbol)==0, NA, as.character(symbol[x]))})
  )

  # Load the gnomAD gene constraints. They have already removed the trailing .version labels for ENST transcripts
  geneConstraints = read.delim(GNOMAD_CONSTRAINT)

  # Search for transcripts/genes in order of ENST, RefSeq, Symbol
  geneConstraints_ind = sapply(1:nrow(query_df), function(x){

    # Check against ENST
    if(!is.na(query_df$ENST_split[x])) {
      x_ind = which(geneConstraints$transcript == query_df$ENST_split[x])
      if(length(x_ind)==1)
        return(x_ind)
    }
    # Check against gene symbol
    if(!is.na(query_df$Symbol[x])) {
      x_ind = which(geneConstraints$gene == query_df$Symbol[x])
      if(length(x_ind)==1)
        return(x_ind)
    }
    return(NA)
  })

  return(geneConstraints[geneConstraints_ind, ])
}


get_gene_chdAnnotation = function(enst=NULL, refseq=NULL, symbol=NULL) {

  query_length = max(length(enst), length(refseq), length(symbol))
  if(!all(c(length(enst), length(refseq), length(symbol)) %in% c(0, query_length)))
    stop("Input query variables have different lengths")

  # Merge input queries into a single data.frame, and remove trailing .version labels for ENST and RefSeq transcripts
  query_df = data.frame(
    "ENST_split" = sapply(1:query_length, function(x){ifelse(length(enst)==0, NA, unlist(strsplit(as.character(enst[x]),"[.]"))[1])}),
    "RefSeq_split" = sapply(1:query_length, function(x){ifelse(length(refseq)==0, NA, unlist(strsplit(as.character(refseq[x]),"[.]"))[1])}),
    "Symbol" = sapply(1:query_length, function(x){ifelse(length(symbol)==0, NA, as.character(symbol[x]))})
  )

  # Load the gene annotation and remove trailing .version labels for ENST and RefSeq transcripts
  library("tools")
  if(file_ext(GENES_CHD_ANNOTATION) == "txt") {
    genesChd = read.delim(GENES_CHD_ANNOTATION)
  } else {
    library("readxl")
    genesChd = read_excel(GENES_CHD_ANNOTATION, 1)
  }
  genesChd$ENST_MANE_Select_split = sapply(as.character(genesChd$ENST_MANE_Select), function(x){unlist(strsplit(x,"[.]"))[1]})
  genesChd$ENST_MANE_PlusClinical_split = sapply(as.character(genesChd$ENST_MANE_PlusClinical), function(x){unlist(strsplit(x,"[.]"))[1]})
  genesChd$RefSeq_MANE_Select_split = sapply(as.character(genesChd$RefSeq_MANE_Select), function(x){unlist(strsplit(x,"[.]"))[1]})
  genesChd$RefSeq_MANE_PlusClinical_split = sapply(as.character(genesChd$RefSeq_MANE_PlusClinical), function(x){unlist(strsplit(x,"[.]"))[1]})
  
  # Search for transcripts/genes in order of ENST, RefSeq, Symbol
  genesChd_ind = sapply(1:nrow(query_df), function(x){

    # Check against ENST MANE Select then MANE Plus Clinical
    if(!is.na(query_df$ENST_split[x])) {
      x_ind = which(genesChd$ENST_MANE_Select_split == query_df$ENST_split[x])
      if(length(x_ind)==1)
        return(x_ind)
      x_ind = which(genesChd$ENST_MANE_PlusClinical_split == query_df$ENST_split[x])
      if(length(x_ind)==1)
        return(x_ind)
    }
    # Check against RefSeq MANE Select then MANE Plus Clinical
    if(!is.na(query_df$RefSeq_split[x])) {
      x_ind = which(genesChd$RefSeq_MANE_Select_split == query_df$RefSeq_split[x])
      if(length(x_ind)==1)
        return(x_ind)
      x_ind = which(genesChd$RefSeq_MANE_PlusClinical_split == query_df$RefSeq_split[x])
      if(length(x_ind)==1)
        return(x_ind)
    }
    # Check against official gene symbol and gene alias
    if(!is.na(query_df$Symbol[x])) {
      x_ind = which(genesChd$Gene == query_df$Symbol[x])
      if(length(x_ind)==1)
        return(x_ind)
      x_ind = which(genesChd$Alias == query_df$Symbol[x])
      if(length(x_ind)==1)
        return(x_ind)
    }
    return(NA)
  })

  return(genesChd[genesChd_ind, ])
}



##### Execute the analyses #####
inputVcfs = paste0(RESULTS_DIR, "/Annotate/annotate_chr", c(1:22,"X"), "_DS02_vep.vcf.gz")
if(!all(file.exists(inputVcfs)))
  stop(paste0("Missing an input vcf at ", RESULTS_DIR, "/Annotate"))


########################### Still to add ###########################
##### Get variant coverage from gnomAD #####
# Keep only rows where VEP symbol matches SpliceAI symbol? What if symbol changed over time?


# Generate tables for various variant annotations, across all chromosomes
variant_kcpra = variant_spliceai = variant_vep = variant_liftover = variant_branchpl = variant_branchpm = variant_lcr = variant_rm = variant_par = variant_sample = NULL
for(vcf in inputVcfs) {
  print(paste("reading in data from:", vcf))
  vcf_variants = get_vcf_variants(vcf)
  if(ncol(vcf_variants) == 0) {
    warning(paste("There are no variants in:", vcf))
    next
  }
  variant_kcpra = rbind(variant_kcpra, vcf_variants)
  variant_spliceai = rbind(variant_spliceai, get_vcf_spliceai(vcf))
  variant_vep = rbind(variant_vep, get_vcf_vepCsq(vcf, mane=MANE_GRCH38))
  variant_liftover = rbind(variant_liftover, get_variant_liftover(vcf=vcf, scratch_dir=paste0(RESULTS_DIR,"/Scratch")))
  variant_branchpl = rbind(variant_branchpl, get_variant_bed(vcf=vcf, bed=BRANCHPOINT_LABRANCHOR_BED_HG38, label="BranchpointL"))
  variant_branchpm = rbind(variant_branchpm, get_variant_bed(vcf=vcf, bed=BRANCHPOINT_MERCER_BED_HG38, label="BranchpointM"))
  variant_lcr = rbind(variant_lcr, get_variant_bed(vcf=vcf, bed=LCR_BED_HG38, label="LCR"))
  variant_rm = rbind(variant_rm, get_variant_bed(vcf=vcf, bed=REPEATMASKER_BED_HG38, label="RepeatMasker"))
  variant_par = rbind(variant_par, get_variant_bed(vcf=vcf, bed=PAR_BED_HG38, label="PAR"))
  variant_sample = rbind(variant_sample, get_variant_altSamples(vcf))
}

# Merge all of the tables
variant_all = Reduce( 
  function(x,y){ merge(x,y,by="Key",all=TRUE) }, 
  list(variant_kcpra, variant_spliceai, variant_vep, variant_liftover, variant_branchpl, variant_branchpm, variant_lcr, variant_rm, variant_par, variant_sample) 
)

# Gene annotations
vaggc = get_gene_gnomadConstraint(enst=variant_all$MANE_ENST, symbol=variant_all$SpliceAI_SYMBOL)
variant_all = cbind(variant_all, vaggc[,c("pLI","oe_lof","oe_lof_lower","oe_lof_upper","oe_mis","oe_mis_lower","oe_mis_upper")])
# vagca = get_gene_chdAnnotation(enst=variant_all$MANE_ENST, refseq=variant_all$Feature, symbol=variant_all$SpliceAI_SYMBOL)
vagca = get_gene_chdAnnotation(enst=variant_all$Feature, refseq=variant_all$MANE_RefSeq, symbol=variant_all$SpliceAI_SYMBOL)
variant_all = cbind(variant_all, vagca[,c("Class_Internal","Class_ClinGen","Subclass","MOI_Internal","MOI_ClinGen","Phenotype1","Phenotype2","OMIM_Phenotypes","OMIM_EvidenceID","OMIM_GeneID","PMID","Remarks","CardiopathySyndromic","ClinGen_Syndrome","ClinGen_AssociatedCardiacDisease","ClinGen_Category")])
variant_all = unique(variant_all)



##### Save the summary to file #####
save(variant_all, file=paste0(RESULTS_DIR, "/Summary/summary_spliceai.rda"))
write.table(variant_all, file=paste0(RESULTS_DIR,"/Summary/summary_spliceai.txt"), sep="\t", quote=F, row.names=F, col.names=T)


