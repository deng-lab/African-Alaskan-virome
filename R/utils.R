library(tidyverse)
library(data.table)
library(ggrepel)

read_CheckV <- function(fpath) {
  dtbl <- fread(fpath) %>%
    mutate(checkv_quality = factor(checkv_quality, levels=c("Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined"))) %>%
    column_to_rownames("contig_id") %>%
    setnames(colnames(.), paste0("checkv_", colnames(.))) %>%
    setnames("checkv_checkv_quality", "checkv_quality") %>%
    rownames_to_column("Contig")
  return(dtbl)
}



read_VirSorter2 <- function(fpath) {
  dtbl <- fread(fpath) %>%
    column_to_rownames("seqname") %>%
    setnames(colnames(.), paste0("virsorter2_", colnames(.))) %>%
    rownames_to_column("Contig")
  return(dtbl)
}



read_vConTACT2 <- function(fpath, assembler_label="_NODE_") {
  # TODO: need to tested using a test dataset
  df_vcontact <- fread(fpath) %>%
    mutate(source=ifelse(str_detect(Genome, assembler_label), "queryseq", "refseq")) %>%
    mutate(cluster_status=ifelse(str_detect(`VC Status`, "Overlap"), "Overlap", `VC Status`)) %>%
    mutate(cluster_status=factor(cluster_status)) %>%
    setnames(colnames(.), paste0("vConTACT_", str_replace_all(colnames(.), " ", "_")))

  vc2_refs <- df_vcontact %>%
    dplyr::filter(str_detect(vConTACT_Genome, assembler_label, negate = TRUE)) %>%
    dplyr::filter(vConTACT_VC != "")

  # Get VC stats
  vcontact_stats <- df_vcontact %>%
    group_by(vConTACT_source, vConTACT_cluster_status) %>%
    dplyr::summarise(seqs_with_VC = sum(vConTACT_VC != ""),
                     seqs_without_VC = sum(vConTACT_VC == "")) %>%
    gather(seq_status, num_seqs, -vConTACT_source, -vConTACT_cluster_status)

  plot_vc_stats <- ggplot(vcontact_stats, aes(y=num_seqs, x=seq_status, fill=vConTACT_source, label=num_seqs)) +
    geom_point(aes(color=vConTACT_source), alpha=0.5, size=3) +
    facet_grid(vConTACT_cluster_status ~ .) +
    # theme(text = element_text(size = 12)) +
    geom_text_repel()

  # Annotate clusters using reference genomes
  vc2_contigs_vclst_anno <- df_vcontact %>%
    # Get distinct cluster ID of contigs in samples
    dplyr::filter(str_detect(vConTACT_Genome, assembler_label)) %>%
    dplyr::filter(vConTACT_VC != "") %>%
    dplyr::select(vConTACT_VC) %>%
    distinct() %>%
    # Annotate cluster with reference taxonomy
    inner_join(vc2_refs, by = "vConTACT_VC")

  # Whether contig clusters were annotated by reference (1) or not (0)
  df_vcontact <- df_vcontact %>%
    mutate(vConTACT_classified=ifelse(vConTACT_VC %in% vc2_contigs_vclst_anno$vConTACT_VC, 1, 0)) %>%
    # only choose assemblies
    dplyr::filter(str_detect(vConTACT_Genome, assembler_label)) %>%
    mutate(vConTACT_VC2=ifelse(vConTACT_VC_Status %in% c("Outlier", "Singleton"), vConTACT_Genome, vConTACT_VC)) %>%
    mutate(vConTACT_VC2=ifelse(str_detect(vConTACT_VC_Status, "Overlap"), vConTACT_VC_Status, vConTACT_VC2)) %>%
    setnames("vConTACT_Genome", "virsorter2_contig_id")

  vc2 <- list("vc_tbl" = df_vcontact,
              "vc_stats" = vcontact_stats,
              "vc_plot" = plot_vc_stats,
              "vc_annotated"= vc2_contigs_vclst_anno,
              "vc_refs" = vc2_refs)
  return(vc2)
}



read_vtaxonomy <- function(fpath) {
  dtbl <- fread(fpath) %>%
    mutate(Contig=str_replace(contig_id, "-cat_[1-6]", "")) %>%
    column_to_rownames("Contig") %>%
    setnames(colnames(.), paste0("vtaxa_", colnames(.))) %>%
    setnames("vtaxa_contig_id", "virsorter2_contig_id") %>%
    mutate(virsorter2_category=as.integer(str_replace(virsorter2_contig_id, ".*-cat_", ""))) %>%
    rownames_to_column("Contig")
  return(dtbl)
}


read_abundance <- function(fpath, anno_viral) {
  # for Phyloseq
  df_abundance <- fread(fpath) %>%
    inner_join(anno_viral[c("Contig", "checkv_contig_length")], by = "Contig") %>%
    arrange(desc(checkv_contig_length), Contig) %>%
    dplyr::select(-checkv_contig_length)
  return(df_abundance)
}


create_viral_abundance <- function(df_abundance, df_covfrac, covfrac_threshold=0.2) {
  df_abundance <- df_abundance %>%
    column_to_rownames("Contig")
  df_filter <- df_covfrac %>%
    column_to_rownames("Contig")

  df_filter[df_filter<covfrac_threshold] <- 0
  df_filter[df_filter>=covfrac_threshold] <- 1
  df_abundance_filtered <- df_abundance*df_filter
  return(df_abundance_filtered)
}


read_lifestyle <- function(fpath) {
  dtbl <- fread(fpath) %>%
    column_to_rownames("V1") %>%
    setnames(colnames(.), paste0("lifestyle_", colnames(.))) %>%
    mutate(lifestyle=ifelse(lifestyle_Temperate>=lifestyle_Virulent, "temperate", "virulent")) %>%
    rownames_to_column("Contig")
  return(dtbl)
}


read_virhost <- function(fpath) {
  # Read only the top hit
  dtbl <- fread(fpath) %>%
    setnames(colnames(.), paste0("vhost_", colnames(.))) %>%
    setnames("vhost_contig_id", "Contig") %>%
    distinct(Contig, .keep_all=T)
  return(dtbl)
}
