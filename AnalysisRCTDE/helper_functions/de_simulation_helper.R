sim_DEGLAM <- function(de_ground_truth, REPLICATES, de_gene, ref, N_samples, nUMI, common_cell_types, 
                       UMI1, UMI_tot, sigma_init, puck, cell_type_info, beta, UMI_vect,
                       other_methods = F,region_orig = NULL, subset_cells = NULL, regularize_expr = F) {
  puck_mod <- puck
  p_res <- matrix(0, nrow = REPLICATES, ncol = 2)
  z_res <- matrix(0, nrow = REPLICATES, ncol = 2)
  e_res <- matrix(0, nrow = REPLICATES, ncol = 2)
  s_res <- matrix(0, nrow = REPLICATES, ncol = 2)
  dec_res <- matrix(0, nrow = REPLICATES, ncol = 2)
  sing_res <- matrix(0, nrow = REPLICATES, ncol = 2)
  bulk_res <- numeric(REPLICATES)
  type1 <- common_cell_types[1]; type2 <- common_cell_types[2]
  ref_counts <- ref@counts[de_gene,]
  for(j in 1:REPLICATES) {
    if(j %% 5 == 0)
      print(j)
    Y <- numeric(N_samples)
    names(Y) <- 1:N_samples
    if(is.null(region_orig))
      region <- sample(c(rep(0, floor(N_samples / 2)), rep(1, ceiling(N_samples/2)))) 
    else
      region <- region_orig
    sample_index_first <- sample(intersect(which(ref@cell_types == type1) , which(nUMI > 1000))) #which cells to include
    sample_index_second <- sample(intersect(which(ref@cell_types == type2) , which(nUMI > 1000))) 
    sample_index_first <- c(sample_index_first, sample_index_first)
    sample_index_second <- c(sample_index_second, sample_index_second)
    for(ind in 1:N_samples) {#N_samples) {
      cell_ind_1 <- sample_index_first[ind]; cell_ind_2 <- sample_index_second[ind]
      if(region[ind] > 0.5) {
        if(regularize_expr) {
          MULT_1 <- sqrt(de_ground_truth[1]); MULT_2 <- sqrt(de_ground_truth[2]);
        } else {
          MULT_1 <- de_ground_truth[1]; MULT_2 <- de_ground_truth[2];
        }
      } else {
        if(regularize_expr) {
          MULT_1 <- 1/sqrt(de_ground_truth[1]); MULT_2 <- 1/sqrt(de_ground_truth[2]);
        } else {
          MULT_1 <- 1; MULT_2 <- 1
        }
      }
      Y[ind] <- rhyper(1,random_round(ref_counts[cell_ind_1]*MULT_1), -ref_counts[cell_ind_1] + ref@nUMI[cell_ind_1], UMI1[ind]) +
        rhyper(1,random_round(ref_counts[cell_ind_2]*MULT_2), -ref_counts[cell_ind_2] + ref@nUMI[cell_ind_2], UMI_tot - UMI1[ind])
    }
    
    X2 <- cbind(1,region)
    X1 <- matrix(0,nrow = N_samples,ncol=0)
    if(is.null(subset_cells))
      sigma_results <- choose_sigma_gene(sigma_init, Y, X1, X2, beta, UMI_vect, 'direct')
    else { 
      my_ind <- (1:length(Y)) %in% sample(1:length(Y), subset_cells)
      sigma_results <- choose_sigma_gene(sigma_init, Y[my_ind], X1[my_ind,], X2[my_ind,], beta[my_ind,], UMI_vect[my_ind], 'direct')
    }
    res <- sigma_results$res
    #print(res$alpha2[2,1])
    
    e_res[j,] <- res$alpha2[2,1:2] #also called alpha_res
    z_res[j, ] <- res$alpha2[2,1:2]/(sqrt(diag(res$I))[c(2,4)])
    s_res[j, ] <- (sqrt(diag(res$I))[c(2,4)])
    p_res[j, ] <- 2*(1 - pnorm(abs(res$alpha2[2,1:2]/(sqrt(diag(res$I))[c(2,4)]))))
    
    if(other_methods) {
      bulk_res[j] <- log(mean(Y[region==1])) - log(mean(Y[region==0]))
      puck_mod@counts[de_gene,] <- Y
      ref_d <- get_decomposed_data_full_doublet(c(de_gene,'Aldoc'), puck_mod, as.matrix(norm_weights), cell_type_info$info)
      dec_de_2 <- log(mean(ref_d@counts[de_gene, which(ref_d@cell_types == 'Purkinje')[region == 1]]*ref_d@nUMI[which(ref_d@cell_types == 'Purkinje')[region == 1]])/mean(ref_d@nUMI[which(ref_d@cell_types == 'Purkinje')[region == 1]])) -
        log(mean(ref_d@counts[de_gene, which(ref_d@cell_types == 'Purkinje')[region == 0]]*ref_d@nUMI[which(ref_d@cell_types == 'Purkinje')[region == 0]])/mean(ref_d@nUMI[which(ref_d@cell_types == 'Purkinje')[region == 0]]))
      dec_de_1 <- log(mean(ref_d@counts[de_gene, which(ref_d@cell_types == 'Bergmann')[region == 1]]*ref_d@nUMI[which(ref_d@cell_types == 'Bergmann')[region == 1]])/mean(ref_d@nUMI[which(ref_d@cell_types == 'Bergmann')[region == 1]])) -
        log(mean(ref_d@counts[de_gene, which(ref_d@cell_types == 'Bergmann')[region == 0]]*ref_d@nUMI[which(ref_d@cell_types == 'Bergmann')[region == 0]])/mean(ref_d@nUMI[which(ref_d@cell_types == 'Bergmann')[region == 0]]))
      dec_res[j,] <- c(dec_de_1, dec_de_2)
      sing_res[j,] <- c(log(mean(Y[beta[,1] > 0.5 & region == 1])) - log(mean(Y[beta[,1] > 0.5 & region == 0])),
                        log(mean(Y[beta[,2] > 0.5 & region == 1])) - log(mean(Y[beta[,2] > 0.5 & region == 0])))
    }
  }
  return(list(e_res = e_res, z_res = z_res, s_res = s_res, p_res = p_res, bulk_de = bulk_res, dec_de = dec_res, sing_de = sing_res))
}

sub_sample_cell <- function(gene_list, raw.data, cell_index, UMI_sample) {
  sub_sample = sample(rep(rownames(raw.data),raw.data[,cell_index]),UMI_sample,replace=FALSE)
  sub_sample = table(sub_sample)[as.character(gene_list)]
  names(sub_sample) = gene_list
  sub_sample[is.na(sub_sample)] <- 0
  return(sub_sample)
}

random_round <- function(x) {
  if(runif(1) < (x %% 1))
    return(ceiling(x))
  return(floor(x))
}

generate_sim_puck <- function(common_cell_types, gene_list, ref, trials = 77) {
  n_cell_types = length(common_cell_types)
  #trials = 77 # 30
  n_conditions = 13
  boundary = ceiling(n_conditions / 2) # DE occurs from conditions 7-13
  N_samples = (n_cell_types * trials * n_conditions * (n_cell_types - 1))/2
  first_UMI = numeric(N_samples); first_type = character(N_samples); second_type = character(N_samples)
  UMI_tot = 1000; UMI_step = UMI_tot / (n_conditions-1)
  UMI_tot = round(UMI_step * (n_conditions-1)); UMI1_vec = round(0:(n_conditions-1)*UMI_step)
  beads = Matrix(0, nrow = length(gene_list), ncol = N_samples)
  rownames(beads) = gene_list; colnames(beads) = 1:N_samples
  firstbeads = Matrix(0, nrow = length(gene_list), ncol = N_samples)
  rownames(firstbeads) = gene_list; colnames(firstbeads) = 1:N_samples
  secbeads = Matrix(0, nrow = length(gene_list), ncol = N_samples)
  rownames(secbeads) = gene_list; colnames(secbeads) = 1:N_samples
  first_index_list <- numeric(N_samples); second_index_list <- numeric(N_samples);
  names(first_index_list) <- 1:N_samples; names(second_index_list) <- 1:N_samples;
  index = 1
  nUMI = ref@nUMI
  for(i in 1:(n_cell_types-1)) {
    print(paste("Progress",i))
    for(j in (i+1):n_cell_types) {
      print(paste("ProgressSecond",j))
      type1 = common_cell_types[i]; type2 = common_cell_types[j]
      for (condition in 1:n_conditions) {
        UMI1 = UMI1_vec[condition]; UMI2 = UMI_tot - UMI1
        for(t in 1:trials) {
          first_UMI[index] = UMI1; first_type[index] = type1; second_type[index] = type2
          firstInd = sample(intersect(which(ref@cell_types == type1) , which(nUMI > 1000)),1)
          secondInd = sample(intersect(which(ref@cell_types == type2) , which(nUMI > 1000)),1)
          first_index_list[index] <- firstInd; second_index_list[index] <- secondInd
          firstbeads[,index] = as.vector(sub_sample_cell(gene_list, ref@counts, firstInd, UMI1))
          secbeads[,index] = as.vector(sub_sample_cell(gene_list, ref@counts, secondInd, UMI2))
          beads[,index] = firstbeads[,index] + secbeads[,index]
          index = index + 1
        }
      }
    }
  }
  UMI_vect <- rep(UMI_tot,N_samples)
  names(UMI_vect) <- colnames(beads)
  puck <- SpatialRNA(NULL, beads, nUMI = UMI_vect, use_fake_coords = T)
  return(puck)
}