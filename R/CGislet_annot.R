#' Annot CpG islets with gene annotation files
#' @param CGislet GRange object of CpG islets
#' @param annot_file path of annot file
#' @param genome Genome sequences, DNAStringSet object 
#' @param annotDb Name of annotate database, org.HS.eg.db for human
#' @param annot_type Format of annot file, default: "gff3"
#' @param tssRegion Range of promoters
#' @return CGislet with CGislet_type
#' @export

CGislet_gene_annot <- function(CGislet,annot_file,genome,annotDb = NULL,
                               annot_type = "gff3",tssRegion = c(-2500,500)){
  
  my_txdb <- txdbmaker::makeTxDbFromGFF(annot_file)
  
  suppressWarnings(annot_chrominfo <- data.frame(seqinfo(my_txdb)))
  
  seq_name <- levels(seqnames(CGislet$CpGislet))
  
  if (any(grepl("^chr|^Chr",rownames(annot_chrominfo)))){
    if (!any(grepl("^chr|^Chr",seq_name))){
      seqlevels(CGislet$CpGislet) <- paste0(substr(rownames(annot_chrominfo)[1],1,3),seq_name)
    }
  } else {
    if (any(grepl("^chr|^Chr",seq_name))){
      seqlevels(CGislet$CpGislet) <- gsub("^M$","MT",gsub("^chr|Chr","",seq_name))
    }
  }
  
  if (!all(seqlevels(CGislet$CpGislet) %in% rownames(annot_chrominfo))){
    stop("\nSeqnames of genome and annot file do not match")
  }
  
  if (is.null(annotDb)){
    islet_Anno <- ChIPseeker::annotatePeak(CGislet$CpGislet, tssRegion=tssRegion, TxDb=my_txdb)
    CGislet$CpGislet <- islet_Anno@anno
    if (annot_type == "gff3" & grepl("gff",gsub(".*/","",annot_file))){
      my_annot <- annot_from_gff3(CGislet,gff3 = annot_file)
      CGislet$CpGislet$SYMBOL <- my_annot[CGislet$CpGislet$geneId,"SYMBOL"]
    }
  } else {
    islet_Anno <- ChIPseeker::annotatePeak(CGislet$CpGislet, tssRegion=tssRegion,
                                           TxDb=my_txdb, annoDb = annotDb)
    CGislet$CpGislet <- islet_Anno@anno
  }
  CGislet$gene_annot <- islet_Anno
  
  seqlevels(CGislet$CpGislet) <- seq_name
  return(CGislet)
  
}


#' Get annotation from annot_list
#' @param annot_list annotation list split from GFF3 annotation column (V9)
#' @param pattern character string
#'

get_annot <- function(annot_list,pattern){
  index <- grep(pattern,annot_list)
  locate <- sapply(annot_list[index],function(x)grep(pattern,x))
  valid <-  annot_list[index]
  result <- gsub(pattern,"",sapply(1:length(valid),
                                   function(x)valid[[x]][locate[x]]))
  return(list(result = result,index = index))
}



#' Get annotation dataframe from gff3 format files
#' @param CGislet CpG islet object
#' @param gff3 gff3 file path of data.frame
#' @param symbol_pattern character;prefix of SYMBOL in annotation file
#' @param ensembl_pattern character;prefix of ENSEMBL names in annotation file
#' @param description logical;if T, return dataframe contain description column
#' @param description_pattern character;prefix of gene descriptions in annotation file
#' @return annotation dataframe includes ENSEMBL, SYMBOL (and descriptions)
#' @export
#'
#'

annot_from_gff3 <- function(CGislet,gff3,symbol_pattern = "Name=",
                            ensembl_pattern = "gene_id=",
                            description = T,
                            description_pattern = "description="){

  if ("character" %in% class(gff3)){
    my_gff3 <- read.table(gff3,sep = "\t",quote = "")
  }

  if ("data.frame" %in% class(gff3)){
    my_gff3 <- gff3
  }

  seq_name <- levels(seqnames(CGislet$CpGislet))

  if (any(grepl("^chr|^Chr",my_gff3$V1))){
    if (!any(grepl("^chr|^Chr",seq_name))){
      seqlevels(CGislet$CpGislet) <- paste0(substr(rownames(annot_chrominfo)[1],1,3),seq_name)
    }
  } else {
    if (any(grepl("^chr|^Chr",seq_name))){
      seqlevels(CGislet$CpGislet) <- gsub("^M$","MT",gsub("^chr|Chr","",seq_name))
    }
  }

  if (!all(seqlevels(CGislet$CpGislet) %in% my_gff3$V1)){
    stop("\nSeqnames of genome and annot file do not match")
  }

  annot <- my_gff3[my_gff3$V3 == "gene",]

  annot <- annot$V9
  annot <- strsplit(annot,split = ";")

  # symbol_index <- grep("Name=",annot)
  # symbol_locate <- sapply(annot[symbol_index],function(x)grep("Name=",x))
  # annot_symbol <- annot[symbol_index]
  # symbol <- gsub("Name=","",sapply(1:length(annot_symbol),
  #                  function(x)annot_symbol[[x]][symbol_locate[x]]))

  ensembl <- get_annot(annot,ensembl_pattern)

  symbol <- get_annot(annot,symbol_pattern)

  result <- data.frame(ENSEMBL = ensembl$result,
                       SYMBOL = "")
  result$SYMBOL[symbol$index] <- symbol$result

  if (description){
    description <- get_annot(annot,description_pattern)
    result$DESCRIPTION <- ""
    result$DESCRIPTION[description$index] <- description$result
  }

  rownames(result) <- result$ENSEMBL

  return(result)
}



#' Get regulatory features annotation of CpG islets
#' @param CG_islet CpG islet object
#' @param regulate_gff3 Regulatory features annotation file; gff3 format
#' @return a list includes regulatory features annotation of CpG islets
#' @export
#'
#'

range_annot <- function(CG_islet,regulate_gff3){

  my_gff3 <- read.table(regulate_gff3,sep = "\t",quote = "")

  seq_name <- levels(seqnames(CG_islet))

  if (any(grepl("^chr|^Chr",my_gff3$V1))){
    if (!any(grepl("^chr|^Chr",seq_name))){
      seqlevels(CG_islet) <- paste0(substr(rownames(annot_chrominfo)[1],1,3),seq_name)
    }
  } else {
    if (any(grepl("^chr|^Chr",seq_name))){
      seqlevels(CG_islet) <- gsub("^M$","MT",gsub("^chr|Chr","",seq_name))
    }
  }

  if (!all(seqlevels(CG_islet) %in% my_gff3$V1)){
    stop("\nSeqnames of genome and annot file do not match")
  }

  my_gff3_list <- split(my_gff3,my_gff3$V3)

  regulate_annot <- rep("",length(CG_islet))

  for (i in 1:length(my_gff3_list)){

    regulate_gff3 <- my_gff3_list[[i]]

    regulate_GRange <- GRanges(seqnames = regulate_gff3$V1,ranges = IRanges(start = regulate_gff3$V4,end = regulate_gff3$V5))

    tmp <- findOverlaps(CG_islet, regulate_GRange)

    tmp_index <- unique(tmp@from)

    regulate_annot[tmp_index] <- paste(names(my_gff3_list)[i],regulate_annot[tmp_index],sep = ";")

  }
  regulate_annot <- strsplit(regulate_annot,";")
  return(regulate_annot)

}

#' Get CGI annotation of CpG islets
#' @param CG_islet CpG islet object
#' @param CGI CGI annotation file, bed format.
#' @return a list includes CGI annotation of CpG islets
#' @export
#'
#'
CGI_annot <- function(CG_islet,CGI){

  my_gff3 <- read.table(CGI,sep = "\t",quote = "")

  # my_gff3_list <- split(my_gff3,my_gff3$V3)

  # regulate_annot <- rep("",length(CG_islet))

  my_seqnames <- unique(my_gff3$V1)

  CG_islet_seqnames <- levels(seqnames(CG_islet))

  if (any(grepl("^chr|^Chr",CG_islet_seqnames))){
    if (!any(grepl("^chr|^Chr",my_seqnames))){
      my_gff3$V1 <- paste0(substr(CG_islet_seqnames[1],1,3),my_gff3$V1)
    }
  } else {
    if (any(grepl("^chr|^Chr",my_gff3$V1))){
      my_gff3$V1 <- gsub("^chr|Chr","",my_gff3$V1)
      my_gff3$V1[which(my_gff3$V1 == "M")] <- "MT"
    }
  }

  regulate_annot <- rep("",length(CG_islet))

  CGI_GRange <- GRanges(seqnames = my_gff3$V1,ranges = IRanges(start = my_gff3$V2,end = my_gff3$V3))

  tmp <- findOverlaps(CG_islet, CGI_GRange)

  tmp_index <- unique(tmp@from)

  regulate_annot[tmp_index] <- paste("CGI",regulate_annot[tmp_index],sep = ";")

  CGI_shores_GRange <- GRanges(seqnames = rep(my_gff3$V1,2),
                               ranges = IRanges(start = c(my_gff3$V2-2000,my_gff3$V3+1),end = c(my_gff3$V2-1,my_gff3$V3+2000)))

  tmp <- findOverlaps(CG_islet, CGI_shores_GRange)

  tmp_index <- unique(tmp@from)

  regulate_annot[tmp_index] <- paste("CGI_shores",regulate_annot[tmp_index],sep = ";")

  CGI_shelves_GRange <- GRanges(seqnames = rep(my_gff3$V1,2),
                               ranges = IRanges(start = c(my_gff3$V2-4000,my_gff3$V3+2001),end = c(my_gff3$V2-2001,my_gff3$V3+4000)))

  tmp <- findOverlaps(CG_islet, CGI_shelves_GRange)

  tmp_index <- unique(tmp@from)

  regulate_annot[tmp_index] <- paste("CGI_shelves",regulate_annot[tmp_index],sep = ";")

  regulate_annot <- strsplit(regulate_annot,";")
  return(regulate_annot)

}

#' Converts genome coordinates from result file of LiftOver
#' @param CGislet GRange object including
#' @param liftover_file character;bed filepath obtained from LiftOver software
#' @return position box ("chr4:100,001-100,001") of converted genome coordinates
#' @export
#'


liftover <- function(CGislet,liftover_file){

  lift <- read.table(liftover_file)

  lift$lift <- paste(paste(lift$V1,lift$V2,sep = ":"),lift$V3,sep = "-")

  lift <- lift[!duplicated(lift$V4),]

  rownames(lift) <- lift$V4

  lift1 <- lift[,"lift"]

  names(lift1) <- rownames(lift)

  tmp <- data.frame(mcols(CGislet))

  tmp$lift_from <- paste(paste(seqnames(CGislet),start(CGislet) + 1,sep = ":"),end(CGislet),sep = "-")

  tmp$lift_to <- lift1[tmp$lift_from]

  return(tmp$lift_to)
}

#' Converts genome coordinates from result file of LiftOver
#' @param CGislet GRange object
#' @param liftover_file character;bed filepath obtained from LiftOver software
#' @return position box ("chr4:100,001-100,001") of converted genome coordinates
#' @export
#'

get_methy_level <- function(CGislet,methy_dat,methy_genome,methy_levels){

  rownames(methy_dat) <- paste(paste(methy_dat$V1,methy_dat$V2,sep = ":"),methy_dat$V3,sep = "-")

  methy <- methy_dat[,methy_levels]

  names(methy) <- rownames(methy_dat)

  methy_genome <- mcols(CGislet)[,methy_genome]

  result <- methy[methy_genome]

  return(result)
}


#' Get tsDMR annotation of CpG islet
#' @param CG_islet gff3 file path of data.frame
#' @param tsDMR tsDMR region GRange list
#' @return a list includes tsDMR annotation of CpG islets
#' @export
#'
#'

tsDMR_annot <- function(CG_islet,tsDMR){

  tsDMR_annot <- rep("",length(CG_islet))

  for (i in 1:length(tsDMR)){

    tsDMR_GRange <- tsDMR[[i]]

    tmp <- findOverlaps(CG_islet, tsDMR_GRange)

    tmp_index <- unique(tmp@from)

    tsDMR_annot[tmp_index] <- paste(names(tsDMR)[i],tsDMR_annot[tmp_index],sep = ";")

  }
  tsDMR_annot <- strsplit(tsDMR_annot,";")
  return(tsDMR_annot)

}
