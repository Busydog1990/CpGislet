#' Get dist of CpG site to nearest another CpG site
#' @param CG_start Start site of CpG sites.
#' @param genome_len Length of genome sequences.
#' @return Dist of CpG site to nearest another CpG site
#' @export

get_CG_dist <- function(CG_start,genome_len){

  dist_right <- c(CG_start[-1] - CG_start[-length(CG_start)],genome_len - CG_start[length(CG_start)]) - 2

  dist_left <- c(CG_start[1],(CG_start[-1] - CG_start[-length(CG_start)] - 2))

  result <- ifelse(dist_right < dist_left,dist_right,dist_left)

  return(result)
}


#' Get dist of CpG site to nearest left CpG site
#' @param CG_start Start site of CpG sites.
#' @param genome_len Length of genome sequences.
#' @return Dist of CpG site to nearest left CpG site
#' @export

get_CG_dist_left <- function(CG_start,genome_len){

  result  <- c(CG_start[1],(CG_start[-1] - CG_start[-length(CG_start)] - 2))

  return(result)
}

#' Get dist of CpG site to nearest right CpG site
#' @param CG_start Start site of CpG sites.
#' @param genome_len Length of genome sequences.
#' @return Dist of CpG site to nearest right CpG site
#' @export

get_CG_dist_right <- function(CG_start,genome_len){

  result  <- c(CG_start[-1] - CG_start[-length(CG_start)],genome_len - CG_start[length(CG_start)]) - 2

  return(result)
}


#' search CpG islet in genome
#' @param CG_stat Information of CpG sites.
#' @param Nmin int. Minimal CpGs dinucleotide in CpG islets. Default:10.
#' @param Dm int. Maximal distance between 2 adjcent CpG islets. Default:75.
#' @param Mer Logical. If true, CpG islets with dist lower than `Merge_len` would be merged.
#' @param Merge_len The CpG islets with a distance no more than M will be merged as a single CpG islet. Default:150.
#' @param genome genome sequences (DNAStringSet).
#' @param flank Flank width (bp) of each CpG site. Default: 500
#' @param write_site_flank_seq Logical. If true, write flank DNA sequence of each CpG site
#' @param write_CpGislet_seq Logical. If true, write CpG islet DNA sequences
#' @return Information of CpG islets and CpG sites.
#' @export

search_CGislet <- function(CG_stat,Nmin,Dm,Mer,Merge_len,
                           getseq = T,genome = genome,flank = 500,
                           write_site_flank_seq = F,write_CpGislet_seq = F){

  my_seqname <- levels(seqnames(CG_stat))

  CG_stat_list <- split(CG_stat,seqnames(CG_stat))

  CpGislet_result <- GRanges()

  CG_stat_result <- GRanges()

  if (getseq){

    CGislet_seq <- DNAStringSet()

    CG_site_flank_seq <- DNAStringSet()

  }

  for (i in 1:length(CG_stat_list)){

  dist_logi <- CG_stat_list[[i]]$dist <= Dm

  dist_str <- paste0(as.numeric(dist_logi),collapse = "")

  if (grepl("^1",dist_str)){

    CpGislet_start <- c(1,gregexpr("01",dist_str)[[1]] + 1)

  } else {

    CpGislet_start <- gregexpr("01",dist_str)[[1]] + 1
  }

  if (grepl("1$",dist_str)){

    CpGislet_end <- c(gregexpr("10",dist_str)[[1]],nchar(dist_str))

  } else {

    CpGislet_end <- gregexpr("10",dist_str)[[1]]

  }

  CG_num <- CpGislet_end - CpGislet_start + 1

  CpGislet_start_choose <- CpGislet_start[CG_num >= Nmin]

  CpGislet_end_choose <- CpGislet_end[CG_num >= Nmin]

  CG_num_choose <- CG_num[CG_num >= Nmin]

  CpGislet_start_pos <- start(CG_stat_list[[i]])[CpGislet_start_choose]

  CpGislet_end_pos <- end(CG_stat_list[[i]])[CpGislet_end_choose]

  if (Mer){

    CpGislet_dist <- CpGislet_start_pos[-1] - CpGislet_end_pos[-length(CpGislet_end_pos)] - 1

    CpGislet_dist_logi <- which(!CpGislet_dist <= Merge_len)

    CpGislet_start_choose <- CpGislet_start_choose[c(1,CpGislet_dist_logi + 1)]

    CpGislet_end_choose <- CpGislet_end_choose[c(CpGislet_dist_logi,length(CpGislet_end_choose))]

    CpGislet_start_pos <- start(CG_stat_list[[i]])[CpGislet_start_choose]

    CpGislet_end_pos <- end(CG_stat_list[[i]])[CpGislet_end_choose]
  }

  CGislet_name <- vector("character",length = length(CG_stat_list[[i]]))

  for (j in 1:length(CpGislet_start_choose)){

    CGislet_name[CpGislet_start_choose[j]:CpGislet_end_choose[j]] <-

      paste0(my_seqname[i],"_",j)
  }


  CpGislet <- GRanges(seqnames = my_seqname[i],ranges = IRanges(start = CpGislet_start_pos,end = CpGislet_end_pos))

  CpGislet$CGislet_name <- paste0(my_seqname[i],"_",1:length(CpGislet))

  CG_stat_list[[i]]$CGislet_name <- CGislet_name

  suppressWarnings(CpGislet_result <- c(CpGislet_result,CpGislet))

  suppressWarnings(CG_stat_result <- c(CG_stat_result,CG_stat_list[[i]]))

  if (getseq){

    seqs <- DNAStringSet(substring(genome[i],CpGislet_start_pos,CpGislet_end_pos))

    CGislet_seq <- c(CGislet_seq,seqs)

    flank_seq <- DNAStringSet(substring(genome[i],start(CG_stat_list[[i]])-flank,end(CG_stat_list[[i]])+flank))

    CG_site_flank_seq <- c(CG_site_flank_seq,flank_seq)
  }
  # print(i)
  }
  if (getseq){

    CpGislet_result$CpG_OE <- get_CpG_OE_ratio(CGislet_seq)
    CpGislet_result$CpA_TpG_OE <- get_CpA_TpG_OE_ratio(CGislet_seq)
    CpGislet_result$CG_density <- vcountPattern("CG",CGislet_seq) / width(CGislet_seq)
    CpGislet_result$GC_content <- get_GC_content(CGislet_seq)
    names(CGislet_seq) <- CpGislet_result$CGislet_name
    CG_stat_result$CpG_OE <- get_CpG_OE_ratio(CG_site_flank_seq)
    CG_stat_result$CpA_TpG_OE <- get_CpA_TpG_OE_ratio(CG_site_flank_seq)
    CG_stat_result$CG_density <- vcountPattern("CG",CG_site_flank_seq) / width(CG_site_flank_seq)
    CG_stat_result$GC_content <- get_GC_content(CG_site_flank_seq)

  }
  if (!getseq){

    result <- list(CpGislet = CpGislet_result,CG_stat = CG_stat_result)

  } else {



    result <- list(CpGislet = CpGislet_result,
                   CG_stat = CG_stat_result)
                   # CpGislet_seq = CGislet_seq,
                   # CG_flank_seq = CG_site_flank_seq)

  }

  if (write_site_flank_seq){

    writeXStringSet(CG_site_flank_seq,file = "CG_site_flank_seq.fa")

  }

  if (write_CpGislet_seq){

    writeXStringSet(CpGislet_seq,file = "CpGislet_seq.fa")

  }

  return(result)

}

#' get_random_seqs from genome
#' @param genome genome sequences (DNAStringSet).
#' @param random_seq_num Number of random sequences
#' @param random_seq_min Minimum length of random sequence.Default: 200
#' @param random_seq_max Maximum length of random sequence.Default: 800
#' @param cutoff Maximum number of non-DNA_BASES.Default: 100
#' @param seeds Random Number Seed.
#' @return Random sequences extracted from imported genome
#' @export
#'
get_random_seqs <- function(genome,random_seq_num = 100,random_seq_min = 200,
                            random_seq_max = 800,cutoff = 100,seeds = 1234){
  bases_all <- Biostrings::letterFrequency(genome,Biostrings::DNA_BASES)
  if (sum(bases_all) < random_seq_min){
    print("Error: Your genome contains less A|T|C|G than minimal length of your random sequences")
    break
  } else {
    set.seed(seeds)
    result_all <- DNAStringSet()
    while (random_seq_num > 0){

      my_random_start <- sample(1:(Biostrings::width(genome)-random_seq_max),size = random_seq_num*2,replace = T)
      my_random_length <- sample(random_seq_min:random_seq_max,size = random_seq_num*2,replace = T)
      my_random_end <- my_random_start + my_random_length - 1
      result <- DNAStringSet(Biostrings::substring(genome,first = my_random_start,last = my_random_end))
      bases <- Biostrings::letterFrequency(result, Biostrings::DNA_BASES)
      result <- result[rowSums(bases) >= (my_random_length - cutoff)]

      if (length(result) >= random_seq_num){
        result_all <- c(result_all,result[1:random_seq_num])
        random_seq_num <- random_seq_num - length(result)
      } else {
        random_seq_num <- random_seq_num - length(result)
        result_all <- c(result_all,result)
      }
    }
    return(result_all)
  }
}

#' Get CG islet candidates
#' @param genome Input genome sequences.
#' @param Nmin int. Minimal CpGs dinucleotide in CpG islets. Default:10
#' @param Dmax int. Maximal distance between 2 adjcent CpG islets. Default:75
#' @param prob numeric. CG site likelihood threshold, CG sites below this threshold will be filtered out. Default:0.6
#' @param Merge Logical. If TRUE, merge CpG islets candidates
#' @param Merge_length int. The CpG islets with a distance no more than M will be merged as a single CpG islet. Default:150
#' @return List includes CpG islets range and CG sites range
#' @export

get_CG_islet <- function(genome,
                         flank_width = 500,
                         Nmin = 10,
                         Dmax = 75,
                         prob = 0.6,
                         Merge = T,
                         Merge_length = 150){

  CG_site <- Biostrings::vmatchPattern("CG",genome)
  genome_len <- width(genome)
  CG_stat <- GenomicRanges::GRanges(seqnames = rep(names(CG_site),sapply(CG_site,length)),
                                    ranges = BiocGenerics::unlist(CG_site))

  CG_start <- lapply(CG_site,start)

  CG_stat$dist <- unlist(lapply(1:length(CG_start),function(x)get_CG_dist(CG_start[[x]],genome_len[x])))

  CGislet <- search_CGislet(CG_stat,Nmin = Nmin,Dm = Dmax,
                            Mer = Merge,Merge_len = Merge_length,
                            getseq = T,genome = genome,flank = flank_width)

  return(CGislet)
}

#' Calculate CG density in CpG range
#'
#' @param site Genomic location of CG site.
#' @param libs Genomic locations of all CG site.
#' @param flank_left Length of left flank of CpG range. Default:500
#' @param flank_right Length of right flank of CpG range. Default:502
#' @return CG density of input sequence.
#' @export

get_CG_density <- function(site,libs,flank_left=500,flank_right=502){
  num <- sum(libs >= site-flank_left & libs <= site+flank_right)
  density <- num / (flank_left + flank_right - 2)
  return(density)
}

#' get_CpA_TpG_OE_ratio
#'
#' @param sequences Input sequence.
#' @return CpA_TpG_OE_ratio of input sequence
#' @export
#' @examples
#' test_seq <- paste0(sample(c("A","T","C","G"),size = 100,replace = T),collapse = "")
#' get_CpA_TpG_OE_ratio(test_seq)

get_CpA_TpG_OE_ratio <- function(sequences){##Calculation of the CpA_TpG_OE_ratio of each sequences
  sequences <- Biostrings::DNAStringSet(sequences)
  bases <- Biostrings::letterFrequency(sequences, DNA_BASES)
  dinucleotide <- Biostrings::dinucleotideFrequency(sequences)
  CA_site <- dinucleotide[,"CA"]
  TG_site <- dinucleotide[,"TG"]
  ratio <- ((CA_site + TG_site)*rowSums(bases))/(bases[,"C"]*bases[,"A"] + bases[,"T"]*bases[,"G"])
  names(ratio) <- NULL
  ratio[!is.finite(ratio)] <- 0
  return(ratio)
}

#' Calculation of the CpG_OE_ratio of each sequences
#' @param sequences Input sequence.
#' @return CpG_OE_ratio of input sequence
#' @export
#' @examples
#' test_seq <- paste0(sample(c("A","T","C","G"),size = 100,replace = T),collapse = "")
#' get_CpG_OE_ratio(test_seq)


get_CpG_OE_ratio <- function(sequences){
  sequences <- Biostrings::DNAStringSet(sequences)
  bases <- Biostrings::letterFrequency(sequences, DNA_BASES)
  dinucleotide <- Biostrings::dinucleotideFrequency(sequences)
  CG_site <- dinucleotide[,"CG"]
  ratio <- (CG_site*rowSums(bases))/(bases[,"C"]*bases[,"G"])
  names(ratio) <- NULL
  ratio[!is.finite(ratio)] <- 0
  return(ratio)
}

#' Calculation of the GC content of each sequences
#' @param sequences Input sequence.
#' @return GC content of input sequence
#' @export
#' @examples
#' test_seq <- paste0(sample(c("A","T","C","G"),size = 100,replace = T),collapse = "")
#' get_GC_content(test_seq)


get_GC_content <- function(sequence){
  sequence <- Biostrings::DNAStringSet(sequence)
  bases <- Biostrings::letterFrequency(sequence,Biostrings::DNA_BASES)
  ratio <- (bases[,"C"] + bases[,"G"])/rowSums(bases)
  names(ratio) <- NULL
  ratio[!is.finite(ratio)] <- 0
  return(ratio)
}

#' Batch Calculation of CpGislet p-values (Improved Bayesian Posterior Probability Computation)
#'
#' @param candidate_seqs DNAStringSet object containing candidate CpGislet sequences
#' @param random_seqs DNAStringSet object containing genomic background sequences
#' @param prior_strength Prior strength parameter controlling the influence of prior distribution (default: 0.5)
#' @param prior_odds Prior odds ratio, P(H1)/P(H0), representing the prior belief in alternative hypothesis (default: 0.5, indicating prior favors H0)
#' @param seed Random seed for reproducibility
#' @param verbose Logical flag controlling detailed output display
#' @return A list containing: 
#' 1) results dataframe with statistical metrics for each candidate sequence, 
#' 2) feature matrices for candidate and background sequences, 
#' 3) Bayesian parameters, and 
#' 4) summary statistics
#'
#' @description
#' This function implements an improved Bayesian approach for identifying CpG islands from candidate DNA sequences.
#' The method computes traditional frequentist p-values based on Mahalanobis distance while providing enhanced
#' Bayesian posterior probability estimates using Zellner's g-prior approach. The algorithm extracts four genomic
#' features (CpG O/E ratio, CpA/TpG O/E ratio, GC content, and CpG density) from both candidate and background
#' sequences, models the null distribution using multivariate normal approximation, and calculates Bayes factors
#' to estimate the posterior probability that each candidate sequence represents a true CpG island (H1).
#' Empirical calibration ensures consistency between frequentist and Bayesian results.
#'
#' @examples
#' \dontrun{
#' # Load necessary libraries
#' library(Biostrings)
#'
#' # Generate example sequences
#' candidate_seqs <- DNAStringSet(c("ACGTACGT", "CGTCGTCG", "ATATATAT"))
#' random_seqs <- DNAStringSet(rep(c("ACGTACGT", "CGTCGTCG", "ATATATAT"), 10))
#'
#' # Run analysis
#' result <- calculate_pvalues_bayesian_batch_v3(
#' candidate_seqs = candidate_seqs,
#' random_seqs = random_seqs,
#' prior_strength = 0.5,
#' prior_odds = 0.5,
#' seed = 123,
#' verbose = TRUE
#' )
#'
#' # View results
#' head(result$results)
#' }
#'
#' @references
#' Zellner, A. (1986). On assessing prior distributions and Bayesian regression analysis with g-prior distributions.
#' In Bayesian inference and decision techniques: Essays in Honor of Bruno de Finetti (pp. 233-243).
#'
#' Benjamin, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach
#' to multiple testing. Journal of the Royal Statistical Society: Series B (Methodological), 57(1), 289-300.
#'
#' @export
calculate_pvalues_bayesian_batch <- function(candidate_seqs, random_seqs,
                                             prior_strength = 0.5,
                                             prior_odds = 0.5,
                                             seed = 1234,
                                             verbose = TRUE) {
    
    # Check input object format
    if (!inherits(candidate_seqs, "DNAStringSet")) {
        stop("candidate_seqs must be a DNAStringSet object")
    }
    if (!inherits(random_seqs, "DNAStringSet")) {
        stop("random_seqs must be a DNAStringSet object")
    }
    
    set.seed(seed)
    
    n_candidates <- length(candidate_seqs)
    
    if (verbose) {
        cat("=== Calculation of CpGislet pvalue ===\n")
        cat("Number of candidate sequences:", n_candidates, "\n")
        cat("Number of random sequences:", length(random_seqs), "\n")
        cat("Prior strength:", prior_strength, "\n")
        cat("Prior probability ratio:", prior_odds, "\n")
    }
    
    # Calculation of eigenvalue
    calculate_features_vectorized <- function(dna_string_set) {
        n_seqs <- length(dna_string_set)
        features <- matrix(NA, nrow = n_seqs, ncol = 4)
        colnames(features) <- c("CpG_OE", "CpA_TpG_OE", "GC_content", "CG_density")
        
        for (i in 1:n_seqs) {
            seq <- as.character(dna_string_set[i])
            seq_no_n <- gsub("N", "", seq)
            
            if (nchar(seq_no_n) == 0) {
                features[i, ] <- rep(NA, 4)
                next
            }
            
            counts <- list(
                A = stringr::str_count(seq_no_n, "A"),
                C = stringr::str_count(seq_no_n, "C"),
                G = stringr::str_count(seq_no_n, "G"),
                T = stringr::str_count(seq_no_n, "T"),
                CG = stringr::str_count(seq_no_n, "CG"),
                CA = stringr::str_count(seq_no_n, "CA"),
                TG = stringr::str_count(seq_no_n, "TG")
            )
            
            L <- nchar(seq_no_n)
            
            # GC content
            GC_content <- (counts$C + counts$G) / L * 100
            
            # CpG O/E ratio
            if (counts$C > 0 && counts$G > 0) {
                CpG_OE <- (counts$CG * L) / (counts$C * counts$G)
            } else {
                CpG_OE <- 0
            }
            
            # CpG density
            CG_density <- counts$CG / L * 100
            
            # CpA|TpG O/E ratio
            denominator <- (counts$C * counts$A + counts$T * counts$G)
			if (is.null(denominator)){
				denominator <- 0
			}
			if (!is.finite(denominator)){
				denominator <- 0
			}
            if (denominator > 0) {
                CpA_TpG_OE <- ((counts$CA + counts$TG) * L) / denominator
            } else {
                CpA_TpG_OE <- 0
            }
            
            features[i, ] <- c(CpG_OE, CpA_TpG_OE, GC_content, CG_density)
        }
        
        return(features)
    }
    

    if (verbose) cat("Calculation of eigenvalue...\n")
    obs_features <- calculate_features_vectorized(candidate_seqs)
    random_features <- calculate_features_vectorized(random_seqs)
    
    na_rows <- rowSums(is.na(random_features)) > 0
    if (any(na_rows)) {
        random_features <- random_features[!na_rows, ]
    }
    
    if (verbose) {
        cat("Random sequence eigenvalue dimension:", dim(random_features), "\n")
    }
    
    # Calculate the mean and covariance of a random sequence
    random_mean <- colMeans(random_features)
    random_cov <- cov(random_features)
    
    # Calculate the posterior probability of each candidate sequence
    n_candidates <- nrow(obs_features)
    results <- data.frame(
        seq_id = 1:n_candidates,
        seq_length = width(candidate_seqs),
        CpG_OE = obs_features[, "CpG_OE"],
        CpA_TpG_OE = obs_features[, "CpA_TpG_OE"],
        GC_content = obs_features[, "GC_content"],
        CG_density = obs_features[, "CG_density"],
        bayesian_pvalue = rep(NA, n_candidates),
        posterior_prob_h1 = rep(0, n_candidates),
        log_bf = rep(NA, n_candidates),
        mahalanobis_dist = rep(NA, n_candidates),
        stringsAsFactors = FALSE
    )
    
    # Calculate the statistical measures for each candidate sequence
    if (verbose) cat("Calculate the statistical measures for each candidate sequence...\n")
    
    for (i in 1:n_candidates) {
        if (verbose && n_candidates > 100 && i %% 100 == 0) {
            cat(paste("  Processing the", i, "/", n_candidates, "candidate sequence\n"))
        }
        
        obs_feature <- obs_features[i, ]

        if (any(is.na(obs_feature))) {
            next
        }
        
        # Calculate Mahalanobis distance
        diff <- obs_feature - random_mean
        mahalanobis_dist <- sqrt(t(diff) %*% solve(random_cov) %*% diff)
        results$mahalanobis_dist[i] <- mahalanobis_dist
        
        # Calculate p-value (based on chi square distribution)
        results$bayesian_pvalue[i] <- pchisq(mahalanobis_dist^2, df = 4, lower.tail = FALSE)
        
        # Posterior probability calculation
        
        # Bayesian factor based on Mahalanobis distance and prior
        D2 <- mahalanobis_dist^2
        
        # Calculate likelihood ratio (H1 relative to H0)
        # H0~N(μ_random, Σ_random)
        # H1~N(μ_candidate, Σ_random),μ_candidate ≠ μ_random
        
        # For a single observation, the maximum likelihood at H1 is at the observation value
        # The likelihood ratio is: L (H1)/L (H0)=exp (0.5 * D2)
        # Use BIC correction
        n_params_h1 <- 4  # 
        n_obs <- 1 
        
        # Bayesian factor approximation (using BIC approximation)
        # BF ≈ exp(0.5 * (BIC_H0 - BIC_H1))
        # BIC = -2 * log(likelihood) + k * log(n)
        log_lik_h0 <- dmvnorm(obs_feature, mean = random_mean, sigma = random_cov, log = TRUE)
        log_lik_h1 <- 0  # Maximizing the likelihood at the observed value to a constant (standardized)
        
        # Use an alternative distribution
        # H1~N(random_mean, τ^2 * Σ_random)
        # Where, the prior of the effect size controlled by τ
        
        # Simplified Bayesian Factor Calculation
        # Using Zellner's g-prior concept
        # BF = (1 + g)^{(-k/2)} * exp( g/(1+g) * D2/2 )
        # Where g is the prior variance ratio parameter
        
        # Bayesian empirical adjustment
        # g = prior_strength * (n_random / n_candidates)
        n_random <- nrow(random_features)
        g <- prior_strength * (n_random / max(n_candidates, 100))
        g <- min(max(g, 0.1), 10)        
        k <- 4 
        
        # Calculate logarithmic Bayesian factor
        log_bf <- (-k/2) * log(1 + g) + (g/(1 + g)) * (D2/2)
        results$log_bf[i] <- log_bf
        
        # Calculate posterior probability
        # P(H1|data) = (BF * prior_odds) / (1 + BF * prior_odds)
        bf <- exp(log_bf)
        posterior_odds <- bf * prior_odds
        posterior_prob_h1 <- posterior_odds / (1 + posterior_odds)

        pval <- results$bayesian_pvalue[i]
        posterior_prob_h1[is.na(posterior_prob_h1)] <- 0
        if (posterior_prob_h1 < 0.1 && pval < 0.05) {
            adjusted_prob <- 1 / (1 + exp(5 * pval - 2))
            posterior_prob_h1 <- max(posterior_prob_h1, adjusted_prob)
        }
        
        results$posterior_prob_h1[i] <- posterior_prob_h1
    }
    
    # Explain Bayesian factors
    interpret_bf <- function(log_bf) {
        bf <- exp(log_bf)
        if (is.na(bf)) return(NA)
        if (bf < 1/3) return("H0")
        else if (bf < 1) return("Weak support for H0")
        else if (bf < 3) return("Weak support for H1")
        else if (bf < 10) return("Medium support for H0")
        else if (bf >= 10) return("Strong support for H1")
    }
    
    results$bf_interpretation <- sapply(results$log_bf, interpret_bf)
    
    # Add significance markers (based on original p-value and posterior probability)
    results$significant_pvalue <- results$bayesian_pvalue < 0.05
    results$significant_prob <- results$posterior_prob_h1 > 0.5  # 
    
    n_sig_prob <- sum(results$significant_prob, na.rm = TRUE)
    n_sig_pval <- sum(results$significant_pvalue, na.rm = TRUE)
    
    # Dynamically adjust the threshold to match the number of p-value filters
    if (n_sig_prob < n_sig_pval * 0.5) {
        prob_threshold <- quantile(results$posterior_prob_h1, 
                                   probs = 1 - (n_sig_pval/n_candidates), 
                                   na.rm = TRUE)
        prob_threshold <- min(max(prob_threshold, 0.1), 0.9)
        results$significant_prob_adjusted <- results$posterior_prob_h1 > prob_threshold
    } else {
        results$significant_prob_adjusted <- results$significant_prob
    }
    
    # Calculate FDR
    results$bayesian_pvalue_adj <- p.adjust(results$bayesian_pvalue, method = "BH")
    
    # Count the number of significant sequences
    n_sig_pval <- sum(results$significant_pvalue, na.rm = TRUE)
	n_sig_pval_adjust <- sum(results$bayesian_pvalue_adj, na.rm = TRUE)
    n_sig_prob <- sum(results$significant_prob, na.rm = TRUE)
	
    if (verbose) {
        cat("\n=== Result Summary ===\n")
        cat("Number of candidates:", n_candidates, "\n")
        cat("Number of Significant candidates (p < 0.05):", n_sig_pval, "(", 
            round(n_sig_pval/n_candidates*100, 1), "%)\n")
		cat("Number of Significant candidates (p.adjust < 0.05):", n_sig_pval, "(", 
            round(n_sig_pval/n_candidates*100, 1), "%)\n")
        cat("Number of Significant candidates (posterior > 0.5):", n_sig_prob, "(", 
            round(n_sig_prob/n_candidates*100, 1), "%)\n")
        results$log_bf
        cat("\nBayesian factor distribution:\n")
        bf_summary <- summary(exp(results$log_bf))
        print(bf_summary)
        
        cat("\nPosterior probability distribution:\n")
        prob_summary <- summary(results$posterior_prob_h1)
        print(prob_summary)
        
        cat("\nCorrelation Analysis:\n")
        cor_test <- cor.test(results$bayesian_pvalue, results$posterior_prob_h1, na.rm = TRUE)
        if (cor_test$p.value > 2.2e-16){
		cat("Pearson correlation coefficient between p-value and posterior probability:", 
            round(cor_test$estimate, 3), " (p =", 
            cor_test$p.value, ")\n")}
		else {
		
		}
    }
    
    # Return results
    return(list(
        results = results,
        obs_features = obs_features,
        random_features = random_features,
        bayesian_params = list(
            prior_strength = prior_strength,
            prior_odds = prior_odds,
            random_mean = random_mean,
            random_cov = random_cov,
            n_random = nrow(random_features)
        ),
        summary = list(
            n_candidates = n_candidates,
            n_sig_pval = n_sig_pval,
            n_sig_prob = n_sig_prob,
            n_sig_pval_adj = n_sig_pval_adjust,
            posterior_prob_summary = summary(results$posterior_prob_h1),
            pvalue_summary = summary(results$bayesian_pvalue)
        )
    ))
}



#' Filter CG islet candidates with pvalue
#' @param CGislet GRange object of CpG islets
#' @param genome Input genome sequences.
#' @param breaks CG_var of random DNA sequences were cut with breaks.
#' @param random_seq_num The number of random DNA sequences.
#' @param CG_var Features which were applied to filter CpG islets.
#' @param pvalueCutoff pvalueCutoff
#' @param qvalueCutoff qvalueCutoff
#' @return List includes filtered CpG islets range and CpG sites range
#' @export

CGislet_filter <- function(CGislet_object,BSgenome,genome,random_seq_rate = 0.0001,
                           CG_var = c("CpG_OE","CpA_TpG_OE","GC_content","CG_density"),
                           pvalueCutoff = 0.05,qvalueCutoff = 0.2){


  CGislet <- CGislet_object$CpGislet

  CG_stat <- CGislet_object$CG_stat

  CGislet_seq <- getSeq(BSgenome,CGislet_object$CpGislet)

  # CG_flank_seq <- CGislet_object$CG_flank_seq

  chrom_lengths <- seqlengths(BSgenome)
  
  randomseq <- do.call(c,lapply(1:length(chrom_lengths),
                              function(x)get_random_seqs(hg38_genome[x],random_seq_min = 200,
                                                         random_seq_max = 1000,
                                                         random_seq_num = round(chrom_lengths[x] * random_seq_rate))))

  
  CGislet$CGislet_pvalue <- exp(rowSums(log(CGislet_stat_result + (1/breaks))))

  CGislet$CGislet_pvalue[is.na(CGislet$CGislet_pvalue)] <- 1

  CGislet$CGislet_qvalue <- p.adjust(CGislet$CGislet_pvalue)

  CGislet_filter <- CGislet[CGislet$CGislet_pvalue <= pvalueCutoff & CGislet$CGislet_qvalue <= qvalueCutoff]

  CGislet_site_filter <- CG_stat[CG_stat$CGislet_name %in% CGislet_filter$CGislet_name]

  CGislet_seq_filter <- CGislet_seq[CGislet_filter$CGislet_name]

  CGislet_filter <- list(CGislet = CGislet_filter,CG_stat = CGislet_site_filter,
                         CGislet_seq = CGislet_seq_filter,random_seq_stat = random_seq_stat)

  return(CGislet_filter)

}


#' Convert GC_content, CpG_OE_ratio, CpA_TpG_OE_ratio, CG_density and dist to nearest CG site
#' of each CG site flank sequences to uniform distribution
#' @param CG_mcol Input parameters of CG sites (flank sequences).
#' @export
#' @return standardized CG sites

standard_CG <- function(CG_mcol){
  f <- function(x){
    test1 <- order(unique(x))
    test2 <- table(x)
    test2 <- test2[order(as.numeric(names(test2)))]
    test3 <- cumsum(test2)
    result <- test3[as.character(x)] / length(x)
    return(result)
  }
  CG_dataframe <- DataFrame(sapply(CG_mcol,function(x)f(x)))
  rownames(CG_dataframe) <- 1:nrow(CG_dataframe)
  #rownames(CG_dataframe) <- rownames(CG_mcol)
  return(CG_dataframe)
}


#' Cluster standardized CG sites based on their CpG_OE_ratio, CpA_TpG_OE_ratio, GC content, CG_density and dist with cmeans
#' @param CG_stat GRanges object with CG stat information.
#' @param CG_var CpG_OE_ratio, CpA_TpG_OE_ratio, GC content, CG_density and dist
#' @param seed Seed of the random number generator.
#' @param centers number of cluster centroid.
#' @param standard Logical. Standardize CG sites or not. Default: TRUE.
#' @param order Logical. Order CG clusters or not. Default: TRUE.
#' @param filter Logical. Filter CG site with max dist. Default: TRUE.
#' @param max_dist Filter CG site with max dist. Default: 150.
#' @return CG levels based on their CpG_OE_ratio, CpA_TpG_OE_ratio, CG_density and dist.
#' @export

get_CG_levels <- function(CG_stat,
                          CG_var = c("CpG_OE","GC_content","CpA_TpG_OE","CG_density","dist"),
                          seed = 100,
                          centers = 2,
                          standard = T,
                          order = T){
  set.seed(seed = seed)
  #CG_stat_mcols <- data.table::data.table(data.frame(mcols(CG_stat)[,CG_var]))

  # CG_stat <- CG_stat(mm9,write_seq = F)
  CG_stat_mcols <- mcols(CG_stat)[,CG_var]

  if (standard){
    standard_CG_stat <- standard_CG(CG_stat_mcols)
    colnames(standard_CG_stat) <- paste0("standard_",colnames(standard_CG_stat))
  } else {
    standard_CG_stat <- CG_stat_mcols
  }


  CG_stat_mcols <- data.table::data.table(data.frame(standard_CG_stat))

  task = as_task_clust(CG_stat_mcols)

  dat <- CG_stat_mcols
  graph = po("learner", lrn("clust.cmeans", centers = centers))
  graph2 = po("learner", lrn("classif.ranger",
                             num.trees = 100,
                             importance = "impurity",
                             predict_type = "prob"))

  learner = as_learner(graph)
  CG_prob <- as.data.table(learner$train(task)$predict(task))
  CG_levels <- CG_prob$partition
  CG_stat_mcols$CG_levels <- CG_levels

  task2 = as_task_classif(CG_stat_mcols,target = "CG_levels")
  learner2 = as_learner(graph2)
  learner2$train(task2)
  prediction = learner2$predict(task2)
  acc <- prediction$score(msr("classif.acc"))
  filter = flt("importance", learner = learner2)
  # filter2 = flt("auc")
  feature_score <- as.data.table(filter$calculate(task2))
  # feature_score2 <- as.data.table(filter2$calculate(task2))
  # CG_standard <- data.frame(CG_standard)
  # my_cluster <- kmeans(CG_standard,centers = centers)
  if (order){
    my_stat <- aggregate(dat,list(CG_levels),mean)
    my_stat$standard_CpA_TpG_OE <- -my_stat$standard_CpA_TpG_OE
    threshold <- apply(my_stat[,-1],2,max)
    cluster_order <- order(rowSums(my_stat[,-1] - rep(apply(my_stat[,-1],2,max),each = centers)),decreasing = T)
    CG_levels <- as.factor(CG_levels)
    levels(CG_levels) <- levels(CG_levels)[cluster_order]
    CG_levels <- as.numeric(as.character(CG_levels))
    CG_levels[dat$standard_CpG_OE >= threshold[1] &
                dat$standard_GC_content >= threshold[2] &
                dat$standard_CpA_TpG_OE <= -threshold[3] &
                dat$standard_CG_density >= threshold[4]] <- 1
    CG_prob$partition <- CG_levels
    CG_prob_order <- c(1,2,cluster_order+2)
    CG_prob <- CG_prob[,..CG_prob_order]
    colnames(CG_prob)[3:(3+centers-1)] <- paste0("prob.",1:centers)
  }
  mcols(CG_stat) <- cbind(mcols(CG_stat),standard_CG_stat,CG_prob)
  result <- list(CG_stat = CG_stat,
                 CG_levels = CG_levels,
                 learners = list(cluster = learner,classif = learner2),
                 feature = list(importance = feature_score,
                                cor = cor(data.frame(CG_stat_mcols)),
                                standard_cor = cor(data.frame(standard_CG_stat))),
                 acc = acc)
  return(result)
}

#' Identification of the structure of CG_islet
#' @param CG_islet CG_site object containing CG_islet and CG_islet_CG_levels column.
#' @param min_core_site Minimum number of core sites per core region.
#' @param min_linker_site Minimum number of linker sites per linker region.
#' @return CG_islet object include structure of CG_islet.
#' @export

CG_islet_structure <- function(CG_site,min_core_site = 5,min_linker_site = 5){

  CG_islet_list <- split(CG_site,CG_site$CGislet_name)[unique(CG_site$CGislet_name)]

  #CG_level_list <- lapply(CG_islet_list,function(x)x$CG_islet_CG_levels)

  CG_type <- list()

  for (i in 1:length(CG_islet_list)){

    # if (!i %% 1000){
    #  print(i)
    # }
    if (sum(CG_islet_list[[i]]$CG_islet == 1) < min_core_site){
      CG_type[[i]] <- rep(2,length(CG_islet_list[[i]]))
    } else if (sum(CG_islet_list[[i]]$CG_islet == 2) < min_linker_site){
      CG_type[[i]] <- rep(1,length(CG_islet_list[[i]]))
    } else {
      CG_type[[i]] <- CG_levels_check(CG_islet_list[[i]]$CG_islet,
                                      min_core_site = min_core_site,min_linker_site = min_linker_site)
      #   }
    }

  }
  CG_site$CG_islet_structure <- as.factor(unlist(CG_type))

  return(CG_site)
}


#' Correct the structure of CG_islet
#' @param CG_levels CG_levels of each CG site (1:core,2:flank)
#' @param min_core_site Minimum number of core sites per core region
#' @param min_linker_site Minimum number of linker sites per linker region

CG_levels_check <- function(CG_levels,min_core_site = 5,min_linker_site = 5){

  CG_levels_seq <- paste0(CG_levels,collapse = "")

  CG_levels_check2_sig <- paste0("12{1,",min_linker_site-1,"}1")

  a = 0

  while (1){

    a = a + 1

    if (a > 100){
      break
    }
    CG_levels_check2 <- gregexpr(CG_levels_check2_sig,CG_levels_seq)[[1]]

    if (CG_levels_check2[1] != -1){

      start2 <- CG_levels_check2

      end2 <- CG_levels_check2+attr(CG_levels_check2,"match.length")-1

      for (i in 1:length(start2)){

        substr(CG_levels_seq,start2[i],end2[i]) <- paste0(rep(1,end2[i]-start2[i]+1),collapse = "")

      }
    } else {
      break
    }
  }

  CG_levels_check1_sig <- paste0("21{1,",min_core_site-1,"}2","|^1{1,",min_core_site-1,"}2","|21{1,",min_core_site-1,"}$")

  CG_levels_check1 <- gregexpr(CG_levels_check1_sig,CG_levels_seq)[[1]]

  b = 0

  while (1){

    b = b + 1

    if (b > 100){
      break
    }

    if (CG_levels_check1[1] != -1){

      start1 <- CG_levels_check1

      end1 <- CG_levels_check1+attr(CG_levels_check1,"match.length")-1

      for (i in 1:length(start1)){

        substr(CG_levels_seq,start1[i],end1[i]) <- paste0(rep(2,end1[i]-start1[i]+1),collapse = "")

      }
    } else {
      break
    }
  }

  # flank_left <- gregexpr("^2+",CG_levels_seq)
  #
  # flank_right <- gregexpr("2+$",CG_levels_seq)

  linker_sig <- paste0("12{",min_linker_site,",}1")

  linker <- gregexpr(linker_sig,CG_levels_seq)[[1]]

  if (linker[1] != -1){

    linker_start <- linker + 1

    linker_end <- linker + attr(linker,"match.length") - 1

    for (i in 1:length(linker_start)){

      substr(CG_levels_seq,linker_start[i],linker_end[i]) <- paste0(rep(0,linker_end[i]-linker_start[i]),collapse = "")

    }

  }

  CG_levels_checked <- strsplit(CG_levels_seq,"")[[1]]

  return(CG_levels_checked)
}

#' get type of CpG islets (core, linker, flank)
#' @param CGislet GRange object of CpG islets
#' @param CG_site GRange object of CpGs in CpG islets
#' @return CGislet with CGislet_type
#' @export

get_CGislet_type <- function(CGislet,CG_site){

  CGislet$CGislet_type <- unlist(lapply(split(CG_site,
                                      CG_site$CGislet_name),
                                function(x)min(as.numeric(as.character(x$CG_islet_structure)))))[CGislet$CGislet_name]

  CGislet$CGislet_type[CGislet$CGislet_type == 0] <- "linker"

  CGislet$CGislet_type[CGislet$CGislet_type == 1] <- "core"

  CGislet$CGislet_type[CGislet$CGislet_type == 2] <- "flank"

  return(CGislet)
}

#' get type of CpG islets (core, linker, flank)
#' @param CGislet GRange object of CpG islets
#' @param annotDb Name of annotate database, org.HS.eg.db for human
#' @param organism organism name in KEGG database.
#' @param GO Logical, if T, do GO enrichment analysis.
#' @param KEGG Logical, if T, do KEGG enrichment analysis.
#' @param GO_pvalueCutoff GO enrichment pvalue cutoff 
#' @param GO_qvalueCutoff GO enrichment qvalue cutoff 
#' @param KEGG_pvalueCutoff KEGG enrichment pvalue cutoff
#' @param KEGG_qvalueCutoff KEGG enrichment qvalue cutoff
#' @return GO/KEGG enrichment analysis
#' @export


CGislet_pathway_enrichment <- function(CGislet,annotDb = NULL,
                                       organism = NULL,
                                       GO = T,
                                       KEGG = T,
                                       GO_pvalueCutoff = 0.05,
                                       GO_qvalueCutoff = 0.1,
                                       KEGG_pvalueCutoff = 0.05,
                                       KEGG_qvalueCutoff = 0.1){

  CGislet_list <- split(CGislet$CGislet,CGislet$CGislet$CGislet_type)

  gene_core <- unique(CGislet_list$core$SYMBOL)

  gene_flank <- unique(CGislet_list$flank$SYMBOL)

  gene_linker <- unique(CGislet_list$linker$SYMBOL)

  gene_all <- unique(CGislet$CGislet$SYMBOL)

  gene_core <- gene_core[nchar(gene_core) > 1]

  gene_flank <- gene_flank[nchar(gene_flank) > 1]

  gene_linker <- gene_linker[nchar(gene_linker) > 1]

  gene_all <- gene_all[nchar(gene_all) > 1]

  gene_core_ONLY <- gene_core[!gene_core %in% c(gene_flank,gene_linker)]

  gene_flank_ONLY <- gene_flank[!gene_flank %in% c(gene_core,gene_linker)]

  gene_linker_ONLY <- gene_linker[!gene_linker %in% c(gene_core,gene_flank)]

  gene_muti <- gene_all[!gene_all %in% c(gene_core_ONLY,gene_flank_ONLY,gene_linker_ONLY)]

  genelist <- list(core = gene_core_ONLY,flank = gene_flank_ONLY,linker = gene_linker_ONLY,
                   muti = gene_muti)

  GO_result <- list()

  KEGG_result <- list()

  for (i in 1:length(genelist)){

    my_symbol <- genelist[[i]]

    my_entrezid <- unique(na.omit(AnnotationDbi::select(annotDb,keys = my_symbol,
                                         keytype = "SYMBOL",columns = c("ENTREZID"))$ENTREZID))

    my_entrezid <- my_entrezid[nchar(my_entrezid) > 1]

    if (GO){

    GO_result[[i]] <- enrichGO(my_entrezid,
                      OrgDb = annotDb, ont='ALL',
                      pAdjustMethod = 'BH',pvalueCutoff = GO_pvalueCutoff,
                      qvalueCutoff = GO_qvalueCutoff,keyType = 'ENTREZID')
    }

    if (KEGG){
    KEGG_result[[i]] <- enrichKEGG(my_entrezid,
                                 organism = organism,
                                 pvalueCutoff = KEGG_pvalueCutoff,
                                 qvalueCutoff = KEGG_qvalueCutoff)
    }
  }
  if (GO){names(GO_result) <- names(genelist)}
  if (KEGG) {names(KEGG_result) <- names(genelist)}

  if (is.null(CGislet$Pathway_annot)){
    Pathway_annot <- list()
  } else {Pathway_annot <- CGislet$Pathway_annot}
  Pathway_annot$Genelist <- genelist
  if (GO){Pathway_annot$GO <- GO_result}
  if (KEGG){Pathway_annot$KEGG <- KEGG_result}

  CGislet$Pathway_annot <- list(GO = GO_result,KEGG = KEGG_result,Genelist = genelist)

  return(CGislet)
}

#' Get CG positions in each CGislet region
#' @param CpG_site CG_site objects contain CG_islet column
#' @return CG_site objects contain number, site_pos and site_pos_type columns

CG_site_pos_detect <- function(CpG_site){

  # tmp <- split(CG_site,CG_site$CGislet_name)
  #
  #
  # CpG_site <- start(CG_site)

  CpG_site_start <- start(CpG_site)

  CG_islet_structure <- CpG_site$CG_islet_structure

  if (all(CG_islet_structure == 2)){

    CG_site_pos <-  (CpG_site_start - CpG_site_start[1])/ (CpG_site_start[length(CpG_site_start)] - CpG_site_start[1])

    CG_site_name <- rep("flank",length(CG_islet_structure))

  } else if (all(CG_islet_structure == 1)){

    CG_site_pos <-  (CpG_site_start - CpG_site_start[1])/ (CpG_site_start[length(CpG_site_start)] - CpG_site_start[1])

    CG_site_name <- rep("core",length(CG_islet_structure))

  } else if (!any(CG_islet_structure == 0)){

    boundary_site <- range(which(CG_islet_structure == 1))

    rep_num <- c((boundary_site[1] - 1),(boundary_site[2] - boundary_site[1] + 1),
                 (length(CG_islet_structure) - boundary_site[2]))

    if (boundary_site[1] != 1){

      CG_site_pos_start <- rep(CpG_site_start[c(1,boundary_site[1],boundary_site[2] + 1)],rep_num)

      CG_site_pos_end <- rep(CpG_site_start[c(boundary_site[1] - 1,boundary_site[2],
                                              length(CG_islet_structure))],rep_num)

    } else {

      CG_site_pos_start <- rep(CpG_site_start[c(boundary_site[1],boundary_site[2] + 1)],rep_num[-1])

      CG_site_pos_end <- rep(CpG_site_start[c(boundary_site[2],length(CG_islet_structure))],rep_num[-1])

    }

    CG_site_pos <- (CpG_site_start - CG_site_pos_start) / (CG_site_pos_end - CG_site_pos_start)

    CG_site_name <- rep(c("flank_left","core","flank_right"),rep_num)

  } else {

    CG_islet_structure_seq <- paste0(CG_islet_structure,collapse = "")

    linker_end <- gregexpr("01",CG_islet_structure_seq)[[1]]

    linker_start <- gregexpr("10",CG_islet_structure_seq)[[1]] + 1

    flank_start <- c()

    flank_end <- c()

    if (grepl("^1",CG_islet_structure_seq)){

      core_start <- c(1,linker_end + 1)

    } else {

      flank_start <- c(flank_start,1)

      core_start <- gregexpr("21|01",CG_islet_structure_seq)[[1]] + 1

      flank_end <- c(flank_end,gregexpr("21",CG_islet_structure_seq)[[1]])

    }

    if (grepl("1$",CG_islet_structure_seq)){

      core_end <- c(linker_start - 1,length(CG_islet_structure))

    } else {

      flank_start <- c(flank_start,gregexpr("12",CG_islet_structure_seq)[[1]] + 1)

      flank_end <- c(flank_end,length(CG_islet_structure))

      core_end <- gregexpr("10|12",CG_islet_structure_seq)[[1]]

    }

    if (length(flank_start) > 0){
      names(flank_start) <- rep("flank",length(flank_start))
      names(flank_end) <- rep("flank",length(flank_end))
    }

    names(core_start) <- names(core_end) <- rep("core",length(core_start))

    names(linker_start) <- names(linker_end) <- rep("linker",length(linker_start))

    struc_start <- sort(c(flank_start,core_start,linker_start))

    struc_end <- sort(c(flank_end,core_end,linker_end))

    CG_site_name <- rep(names(struc_start),c(struc_end - struc_start + 1))

    linker_range <- range(which(CG_islet_structure == 0))

    CG_site_name[1:(linker_range[1]-1)] <-
      paste0(CG_site_name[1:(linker_range[1]-1)],"_left")

    CG_site_name[(linker_range[2] + 1):length(CG_site_name)] <-
      paste0(CG_site_name[(linker_range[2] + 1):length(CG_site_name)],"_right")

    CG_site_pos_start <- rep(CpG_site_start[struc_start],c(struc_end - struc_start + 1))

    CG_site_pos_end <- rep(CpG_site_start[struc_end],c(struc_end - struc_start + 1))

    CG_site_pos <- (CpG_site_start - CG_site_pos_start) / (CG_site_pos_end - CG_site_pos_start)


  }

  names(CG_site_pos) <- CG_site_name

  return(CG_site_pos)
}

#' Get CpG positions in each CGislet region
#' @param CGislet CGislet objects contain CGislet_type column
#' @param CpG_site CG_site objects contain CG_islet column
#' @return CG_site with CpG positions
#' @export

CG_site_pos <- function(CGislet,CG_site){

  CG_site_list <- split(CG_site,CG_site$CGislet_name)[CGislet$CGislet_name]

  CG_site_pos <- lapply(CG_site_list,CG_site_pos_detect)

  CG_site$CG_site_pos <- unlist(CG_site_pos)

  CG_site$type <- gsub("^.*\\.","",names(unlist(CG_site_pos)))

  return(CG_site)

}

#' Detection of flank CpG site positions of each CpG islet
#' @param CGislet CGislet objects contain CGislet_type column
#' @param CpG_site CG_site objects contain CG_islet column
#' @param flank Length of left and right flank of each CpG islet.Default: 2000L.
#' @param filter Logical. If TRUE, CpG sites in CpG islet were NOT considered as flank CpG sites.
#' @return CG_site with flank CpG site positions
#' @export

get_CGislet_flank_site_pos <- function(CGislet,CG_stat,flank = 2000L,filter = T){

  CGislet_left <- GRanges(seqnames = seqnames(CGislet),ranges = IRanges(start = start(CGislet) - flank - 1,end = start(CGislet) - 1))

  CGislet_right <- GRanges(seqnames = seqnames(CGislet),ranges = IRanges(start = end(CGislet) + 1,end = end(CGislet) + flank + 1))

  start_CGislet_left <- start(CGislet)- flank - 1

  names(start_CGislet_left) <- CGislet$CGislet_name

  start_CGislet_right <- end(CGislet) + 1

  names(start_CGislet_right) <- CGislet$CGislet_name

  left_flank <- findOverlaps(CG_stat,CGislet_left)

  right_flank <- findOverlaps(CG_stat,CGislet_right)

  CG_stat$left_flank <- ""

  CG_stat$left_flank[from(left_flank)] <- CGislet$CGislet_name[to(left_flank)]

  CG_stat$right_flank <- ""

  CG_stat$right_flank[from(right_flank)] <- CGislet$CGislet_name[to(right_flank)]

  if (filter){

    CG_stat$left_flank[nchar(CG_stat$CGislet_name) > 0] <- ""

    CG_stat$right_flank[nchar(CG_stat$CGislet_name) > 0] <- ""

  }

  CG_stat$left_flank_pos <- as.numeric(((start(CG_stat) - start_CGislet_left[CG_stat$left_flank]) + 1) / (flank + 1))

  CG_stat$right_flank_pos <- as.numeric(((start(CG_stat) - start_CGislet_right[CG_stat$right_flank]) + 1)  / (flank + 1))

  return(CG_stat)

}


#' Extract sequences of CG islets.
#' @param CG_islet CG_islet object.
#' @param CG_site CG_site object.
#' @param genome Input genome sequences.
#' @param range If TRUE. Also return a genome range of choosed region. Valid if region in c("core","linker","flank")
#' @param write_seq if True,write CG islets sequences.Default:FALSE.
#' @param file_name filenames of CG islets sequences.Default:"CG_islet.fa".
#' @param region Obtain the sequence of the specified area in CGislets. One of "all","core","linker" or "flank"
#' @return CG islets ranges with their CpG_OE_ratio CpA_TpG_OE_ratio CG_density and GC_content
#' @export
#'
#'
get_CGislet_seq <- function(CG_islet,CG_site,genome,range = T,write_seq = FALSE,
                            file_name = "CG_islet.fa",region = "all"){

  if (region == "all"){

    CG_islet_list <- split(CG_islet,seqnames(CG_islet))

    CG_islet_seqname <- c()

    CG_islet_seq <- Biostrings::DNAStringSet()

    for (i in 1:length(CG_islet_list)){
      #print(i)
      tmp <- Biostrings::DNAStringSet(substring(genome[names(CG_islet_list)[i]],
                                                start(CG_islet_list[[i]]),
                                                end(CG_islet_list[[i]])))
      CG_islet_seq <- c(CG_islet_seq,tmp)

      CG_islet_seqname <- c(CG_islet_seqname,CG_islet_list[[i]]$CGislet_name)
    }

    names(CG_islet_seq) <- CG_islet_seqname

    if (write_seq){

      Biostrings::writeXStringSet(CG_islet_seq,filepath = file_name)

    }
    return(CG_islet_seq)

  } else if (region == "core"){

    CG_site_list <- split(CG_site,CG_site$CGislet_name)

    CGislet_core <- CG_islet$CGislet_name[CG_islet$CGislet_type == "core"]

    CG_site_core <- CG_site[(CG_site$CG_islet %in% CGislet_core) & CG_site$CG_islet_structure == 1]

    CG_site_core_list <- split(CG_site_core,CG_site_core$CG_islet)

    CG_site_core_start <- lapply(CG_site_core_list,function(x)start(x[1]))

    CG_site_core_end <- lapply(CG_site_core_list,function(x)end(x[length(x)]))

    CG_site_core_dataframe <- data.frame(start = unlist(CG_site_core_start),
                                         end = unlist(CG_site_core_end),
                                         CG_islet = paste0(names(CG_site_core_list),"_core"),
                                         chr = gsub("_.*","",names(CG_site_core_list)))

    CGislet_linker <- CG_islet$CGislet_name[CG_islet$CGislet_type == "linker"]

    CG_site_linker <- CG_site[(CG_site$CG_islet %in% CGislet_linker) & CG_site$CG_islet_structure %in% c(1,3)]

    CG_site_linker_list <- split(CG_site_linker,CG_site_linker$CG_islet)


    for (i in 1:length(CG_site_linker_list)){

      structures <- paste0(CG_site_linker_list[[i]]$CG_islet_structure,collapse = "")
      core_start <- c(1,gregexpr("31",structures)[[1]] + 1)
      core_end <- c(gregexpr("13",structures)[[1]],nchar(structures))
      tmp_data.frame <- data.frame(start = start(CG_site_linker_list[[i]][core_start]),
                                   end = end(CG_site_linker_list[[i]][core_end]),
                                   CG_islet = paste0(names(CG_site_linker_list)[i],"_core"),
                                   chr = gsub("_.*","",names(CG_site_linker_list)[i]))
      tmp_data.frame$CG_islet <- paste0(tmp_data.frame$CG_islet,1:nrow(tmp_data.frame))

      CG_site_core_dataframe <- rbind(CG_site_core_dataframe,tmp_data.frame)
    }

    CG_islet_list <- split(CG_site_core_dataframe,CG_site_core_dataframe$chr)


    CG_islet_seqname <- c()

    CG_islet_seq <- Biostrings::DNAStringSet()

    for (i in 1:length(CG_islet_list)){
      #print(i)
      tmp <- Biostrings::DNAStringSet(substring(genome[names(CG_islet_list)[i]],
                                                CG_islet_list[[i]]$start,
                                                CG_islet_list[[i]]$end))
      CG_islet_seq <- c(CG_islet_seq,tmp)

      CG_islet_seqname <- c(CG_islet_seqname,CG_islet_list[[i]]$CG_islet)
    }

    names(CG_islet_seq) <- CG_islet_seqname

    if (write_seq){

      Biostrings::writeXStringSet(CG_islet_seq,filepath = file_name)

    }

    if (range){

      CG_site_core_range <- GRanges(seqnames = CG_site_core_dataframe$chr,
                                    ranges = IRanges(start = CG_site_core_dataframe$start,
                                                     end = CG_site_core_dataframe$end),
                                    CG_islet = CG_site_core_dataframe$CG_islet)
      result <- list(Range = CG_site_core_range,seq = CG_islet_seq)

      return(result)
    }
    return(CG_islet_seq)

  } else if (region == "flank"){

    CGislet_core <- CG_islet$CGislet_name[CG_islet$CGislet_type != "flank"]

    CG_site_core <- CG_site[CG_site$CG_islet %in% CGislet_core]

    CG_site_list <- split(CG_site_core,CG_site_core$CG_islet)

    structures <- unlist(lapply(CG_site_list,function(x)paste0(x$CG_islet_structure,collapse = "")))

    right_end <- nchar(structures)

    flank_left_end <- unlist(lapply(gregexpr("^2+",structures),function(x)attr(x,"match.length") + 1));names(flank_left_end) <- names(CG_site_list)

    flank_left_start <- ifelse(flank_left_end == 0,NA,1)

    flank_right <- unlist(lapply(gregexpr("2+$",structures),function(x)attr(x,"match.length")));names(flank_right) <- names(CG_site_list)

    flank_right_start <- ifelse(flank_right == -1,NA,right_end-flank_right)

    flank_right_end <- ifelse(flank_right == -1,NA,right_end)

    tmp_data.frame_left <- na.omit(data.frame(start = flank_left_start,
                                              end = flank_left_end,
                                              CG_islet = paste0(names(CG_site_list),"_flank_left"),
                                              chr = gsub("_.*","",names(CG_site_list))))
    CG_site_list_left <- CG_site_list[rownames(tmp_data.frame_left)]



    tmp_data.frame_right <- na.omit(data.frame(start = flank_right_start,
                                               end = flank_right_end,
                                               CG_islet = paste0(names(CG_site_list),"_flank_right"),
                                               chr = gsub("_.*","",names(CG_site_list)),
                                               row.names = names(CG_site_list)))

    CG_site_list_right <- CG_site_list[rownames(tmp_data.frame_right)]

    for (i in 1:nrow(tmp_data.frame_left)){

      #print(i)
      tmp_data.frame_left$start[i] <- Biostrings::start(CG_site_list_left[[i]][tmp_data.frame_left$start[i]])

      tmp_data.frame_left$end[i] <- end(CG_site_list_left[[i]][tmp_data.frame_left$end[i]])

    }

    for (i in 1:nrow(tmp_data.frame_right)){

      tmp_data.frame_right$start[i] <- start(CG_site_list_right[[i]][tmp_data.frame_right$start[i]])

      tmp_data.frame_right$end[i] <- end(CG_site_list_right[[i]][tmp_data.frame_right$end[i]])

    }

    CG_site_flank_dataframe <- rbind(tmp_data.frame_left,tmp_data.frame_right)

    CGislet_flank <- CG_islet[CG_islet$CGislet_type == "all_flank"]

    CG_site_flank <- data.frame(start = start(CGislet_flank),
                                end = end(CGislet_flank),
                                CG_islet = paste0(CGislet_flank$CGislet_name,"_flank"),
                                chr = gsub("_.*","",CGislet_flank$CGislet_name))
    CG_site_flank_dataframe <- rbind(CG_site_flank_dataframe,CG_site_flank)

    CG_islet_list <- split(CG_site_flank_dataframe,CG_site_flank_dataframe$chr)

    CG_islet_seqname <- c()

    CG_islet_seq <- Biostrings::DNAStringSet()

    for (i in 1:length(CG_islet_list)){
      #print(i)
      tmp <- Biostrings::DNAStringSet(substring(genome[names(CG_islet_list)[i]],
                                                CG_islet_list[[i]]$start,
                                                CG_islet_list[[i]]$end))
      CG_islet_seq <- c(CG_islet_seq,tmp)

      CG_islet_seqname <- c(CG_islet_seqname,CG_islet_list[[i]]$CG_islet)
    }

    names(CG_islet_seq) <- CG_islet_seqname

    if (write_seq){

      Biostrings::writeXStringSet(CG_islet_seq,filepath = file_name)

    }

    if (range){

      CG_site_flank_range <- GRanges(seqnames = CG_site_flank_dataframe$chr,
                                     ranges = IRanges(start = CG_site_flank_dataframe$start,
                                                      end = CG_site_flank_dataframe$end),
                                     CG_islet = CG_site_flank_dataframe$CG_islet)
      result <- list(Range = CG_site_flank_range,seq = CG_islet_seq)

      return(result)
    }

    return(CG_islet_seq)

  } else if (region == "linker"){

    CGislet_linker <- CG_islet$CGislet_name[CG_islet$CGislet_type == "linker"]

    CG_site_linker <- CG_site[(CG_site$CG_islet %in% CGislet_linker) & CG_site$CG_islet_structure %in% c(1,3)]

    CG_site_linker_list <- split(CG_site_linker,CG_site_linker$CG_islet)

    CG_site_linker_dataframe <- data.frame()

    for (i in 1:length(CG_site_linker_list)){

      structures <- paste0(CG_site_linker_list[[i]]$CG_islet_structure,collapse = "")
      linker_start <- gregexpr("13",structures)[[1]]
      linker_end <- gregexpr("31",structures)[[1]] + 1
      tmp_data.frame <- data.frame(start = start(CG_site_linker_list[[i]][linker_start]),
                                   end = end(CG_site_linker_list[[i]][linker_end]),
                                   CG_islet = paste0(names(CG_site_linker_list)[i],"_linker"),
                                   chr = gsub("_.*","",names(CG_site_linker_list)[i]))
      tmp_data.frame$CG_islet <- paste0(tmp_data.frame$CG_islet,1:nrow(tmp_data.frame))

      CG_site_linker_dataframe <- rbind(CG_site_linker_dataframe,tmp_data.frame)
    }

    CG_islet_list <- split(CG_site_linker_dataframe,CG_site_linker_dataframe$chr)

    CG_islet_seqname <- c()

    CG_islet_seq <- Biostrings::DNAStringSet()

    for (i in 1:length(CG_islet_list)){
      #print(i)
      tmp <- Biostrings::DNAStringSet(substring(genome[names(CG_islet_list)[i]],
                                                CG_islet_list[[i]]$start,
                                                CG_islet_list[[i]]$end))
      CG_islet_seq <- c(CG_islet_seq,tmp)

      CG_islet_seqname <- c(CG_islet_seqname,CG_islet_list[[i]]$CG_islet)
    }

    names(CG_islet_seq) <- CG_islet_seqname

    if (write_seq){

      Biostrings::writeXStringSet(CG_islet_seq,filepath = file_name)

    }

    if (range){

      CG_site_linker_range <- GRanges(seqnames = CG_site_linker_dataframe$chr,
                                      ranges = IRanges(start = CG_site_linker_dataframe$start,
                                                       end = CG_site_linker_dataframe$end),
                                      CG_islet = CG_site_linker_dataframe$CG_islet)
      result <- list(Range = CG_site_linker_range,seq = CG_islet_seq)

      return(result)
    }
    return(CG_islet_seq)



  } else {

    stop("invalid CGislet type")

  }
}

#' Permutation test of features which were applied to filter CpG islets.
#' @param dat CG_islet list.
#' @param n The number of permutations in permutation testing.
#' @param vars Features which were applied to filter CpG islets.
#' @return pvalue of permutation_test
#' @export
#'

permutation_test <- function(dat,n = 1000,vars = "CpG_OE"){
  dat <- split(dat,dat$CGislet_type)
  result <- matrix(0,length(dat),length(dat))
  for (i in 1:length(dat)){
    for (j in 1:length(dat)){
      if (i < j){
        pvalue <- wilcox.test(dat[[i]][,vars],dat[[j]][,vars])$p.value
        pvalue_per <- c()
        dat_all <- rbind(dat[[i]],dat[[j]])
        for (k in 1:n){
          tmp_sample <- sample(1:nrow(dat_all),size = nrow(dat[[i]]), replace = F)
          pvalue_tmp <- wilcox.test(dat_all[tmp_sample,vars],dat_all[-tmp_sample,vars])$p.value
          pvalue_per <- c(pvalue_per,pvalue_tmp)
        }
        result[i,j] <- sum(pvalue > pvalue_per) / n
      }
    }
  }
  return(result)
}

#' get pvalue with wilcox.test and permutation_test compared two group of features which were applied to filter CpG islets.
#' @param dat CG_islet list.
#' @param bootstrap Bootstrap sampling frequency.Default: 100.
#' @param permutation The number of permutations in permutation testing.
#' @param vars Features which were applied to filter CpG islets.
#' @param seed Seed of the random number generator.
#' @return pvalue of two group of features
#' @export
#'
get_pvalue <- function(dat,bootstrap = 100L,permutation = 1000L,vars = "CpG_OE",seed = 1234){

  set.seed(1234)

  dat_list <- split(dat,f = dat$CGislet_type)

  pvalue_list <- vector("list", bootstrap)

  for (i in 1:bootstrap){

    dat_bootstrap <- do.call(rbind,lapply(dat_list,
                                          function(x)x[sample(1:nrow(x),size = bootstrap,replace = T),]))

    pvalue_list[[i]] <- permutation_test(dat_bootstrap,n = permutation,vars = vars)

  }

  pvalue <- Reduce("+",pvalue_list)
  pvalue <- pvalue / bootstrap


  colnames(pvalue) <- rownames(pvalue) <- names(dat_list)

  return(pvalue)

}


