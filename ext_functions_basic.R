#................................................
# Procedure: Functions for basic TCI algorithm
# Project: TCI
# User: LXAY
# Date: Apr-2021
#................................................


  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Function to compute scores with the Basic algorithm (tci, probability, confidence)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  
  basic_algorithm = function(  tdatasets = NULL, # Table of datasets
                                genes.query = NULL, # List of query genes
                                ldatasets = NULL, # List datasets (matrices)
                                tscore.attr = NULL, # Attributes for computing the score
                                plot.file.basic = NULL, # Path to plot file
                                n.permutations = 10000, # Number of permutations
                                save.data.objects.folder = NA, # Folder to save data objects (e.g., raw and transformed datasets)
                                permute.option = c("rows and columns", "columns", "rows", "rows per layer", "rows per column") # Option to permute rows, columns or both in evidence.trans.matrix
          ){
    
    
    cat("-------------------\n Running basic algorithm for scoring genes. \n------------------- \n")
    
    
    
    
    # Select evidence layers 
    levidence.layers = tlayers$evidence_layer [ tlayers$evidence_layer_weight != "no" ]
    
    
    # If the user selected "confidence score" to be computed
    if( tscore.attr$compute [ tscore.attr$score_attribute == "confidence" ] == "yes" ){
      
      # Initialise table of confidence scores per layer
      tconfidence.score = NULL
     
      
      # Call function to compute confidence value based on data completeness
      cat("-------------------\n Call function to compute confidence value based on data completeness. \n------------------- \n")
      # Loop across evidence layers
      for(i in 1:length(levidence.layers)){
        
        # If there are not datasets for this layer shows an error
        if ( sum( tdatasets$evidence_layer == levidence.layers[i] ) == 0 ){
          stop( paste0("There are not datasets uploaded for this evidence layer:", levidence.layers[i])  )
        }
        
        
        # Call function for all datasets  
        tconfidence.score = cbind(tconfidence.score,  
                                  unlist( get( "gene.confidence.score", confidence_score( ldatasets = ldatasets[which( tdatasets$evidence_layer == levidence.layers[i] ) ] , genes.query = genes.query,
                                                                                          main = paste("layer: ", levidence.layers[i], sep = "" ))
                                                                                          ) ) )
                                                                
        
        
      }
      
      head(tconfidence.score)
      
      # Compute the total confidence score given the set of weights
      cat("-------------------\n Compute the total confidence score given the set of weights for the layers. \n------------------- \n")
      
      tconfidence.score = cbind( genes.query,
                                 tconfidence.score,
                                 round( tconfidence.score %*% as.numeric( tlayers$evidence_layer_weight [ tlayers$evidence_layer_weight != "no" ]), 3) ) 
      
      
      
      # Add colnames to table of confidence values
      colnames(tconfidence.score) <- c("ensembl_gene_id", paste("confidence_score_", levidence.layers, sep=""), "confidence_score_total")
      rownames(tconfidence.score) <- NULL
      head(tconfidence.score)
    
      
    }else{
      
      tconfidence.score = NA
      
    } # End if the user selected "confidence score" to be computed
    
    
    # The basic score do not need data imputation step
    
    
    # Transform the data to the interval 0,1
    cat("-------------------\n Transform the data and scale to the interval 0,1. \n------------------- \n")
    
  
    res.basic.transform.data = basic_transform_data(ldatasets = ldatasets,
                                                    tdatasets = tdatasets)
    
    # Extract results from data transformation
    ltrans.datasets <- res.basic.transform.data$ltrans.datasets       # lapply(ltrans.datasets, nrow)     # lapply(ldatasets, nrow)
    
    # Assign names to datasets
    names(ltrans.datasets) <- names(ldatasets)
    
    # Put all transformed datasets in one single matrix
    res.merge.datasets = merge_datasets(ldatasets = ltrans.datasets [ which( !is.na( tdatasets$weights) )],
                                        genes.query = genes.query,
                                        summarise = tdatasets$curated_dataset_columns_agregation [ which( !is.na( tdatasets$weights) )]) 
    
    # Extract resuls from data merging
    evidence.trans.matrix = res.merge.datasets$evidence.matrix
    
    # Transform columns to numeric
    cat(" ++ Transforming columns to numeric \n")
    row.names.matrix <- rownames(evidence.trans.matrix)
    evidence.trans.matrix <- apply(evidence.trans.matrix, MARGIN = 2, FUN = as.numeric)
    rownames(evidence.trans.matrix) <- row.names.matrix
    evidence.trans.matrix = as.data.frame(evidence.trans.matrix, stringsAsFactors = F)
    print( summary(evidence.trans.matrix) )
    
    
    # Print progress
    cat("++ Gene scores per dataset: \n")
    
    print(head(evidence.trans.matrix ))
    
    
    # print(ncol( evidence.trans.matrix ) );   print(length(tdatasets$weights))
    if( ncol( evidence.trans.matrix ) != sum( !is.na( tdatasets$weights) )) {
      stop( "In the merged matrix of transformed datasets: More columns than weights.")
    }
    
    
    
    # Initialise table of tci per layer
    ttci.score = NULL
    
    
    # Check if the user selected "P(S>s) for tci" to be computed
    p.value.permutation = tscore.attr$compute [ tscore.attr$score_attribute == "randomisation" ] == "yes"
    
    
    # Initialise table of P(S>s)s of tci per layer 
    # Number of permutations where the gene achieved a score higher than the one observed
    
    if( p.value.permutation == TRUE ){ 
      ttci.probability = NULL 
    }else{
      ttci.probability = NA
    }
    
    
    # Compute tci evidence score per evidence layer
    cat("-------------------\n Compute tci evidence score per evidence layer. \n------------------- \n")
  
    par(mfrow = c(2,2))
    
    for(i in 1:length(levidence.layers)){
      
      # Call tci function for the datasets in this layer
      cat(" ++ Calling tci function for the datasets in layer:", levidence.layers[i], "\n")
      
      # Multiply each dataset for a given weight
      dim(evidence.trans.matrix)
      
      # Sum the scores across data sets per layer and compute the layer score (0-1) and total score.
      res.basic.tci.score = basic_tci_score( evidence.trans.matrix = evidence.trans.matrix[  , tdatasets$evidence_layer == levidence.layers[i], drop = FALSE  ] ,
                                             basic.weights = tdatasets$weights [ tdatasets$evidence_layer == levidence.layers[i] ],
                                             main = paste("TCI, layer: ", levidence.layers[i], sep = "" ),
                                             replace.NA.0 = TRUE, 
                                             p.value.permutation = FALSE, 
                                             n.permutations = NULL )
      
      # Add score vector to table of tci                      
      ttci.score = cbind(  ttci.score,  
                           res.basic.tci.score$gene.tci.score )
      
     
      
      
      
    } # End loop of tci evidence score per evidence layer
    
    
    # Compute the total tci
    cat("-------------------\n Compute the total tci. \n------------------- \n")
    res.basic.tci.score = basic_tci_score( evidence.trans.matrix = ttci.score,
                                           basic.weights = as.numeric( tlayers$evidence_layer_weight [ tlayers$evidence_layer_weight != "no" ], na.rm = T), 
                                           main = "Total TCI", 
                                           replace.NA.0 = TRUE,  
                                           p.value.permutation = FALSE, # The probabilities are computed below using the whole evidence matrix
                                           n.permutations = NULL )
    
    
    
    # Add the total tci to the table of tci
    ttci.score = cbind(rownames(evidence.trans.matrix), 
                       ttci.score, 
                       res.basic.tci.score$gene.tci.score)
    
    
    # Add colnames to the table and transform to data frame
    colnames(ttci.score) <- c("ensembl_gene_id", paste("tci_", levidence.layers, sep = ""), "tci_total")
    ttci.score = data.frame(ttci.score, stringsAsFactors = F)
    print( head(ttci.score) )  # sum(is.na(ttci.score$total_tci))
    rownames(ttci.score) <- rownames(evidence.trans.matrix)
    
    # Prepare table of p.values of tci
    if( p.value.permutation == TRUE ){ 
      
      # Initialise table of randomness score
      ttci.probability = rownames(ttci.score)
      
      # Compute the P(S>s) for each layer
      for( i in 1:nrow(tlayers)){
        
        cat(" ++ Computing P(S>s) for layer", tlayers$evidence_layer[i], ". \n")
        
        
        # Use the evidence matrix to permute the scores (only columns of the layer)
        pos.layer = colnames(evidence.trans.matrix) %in% tdatasets$curated_dataset_id[tdatasets$evidence_layer == tlayers$evidence_layer[i]] 
        res.basic.tci.score = basic_tci_score( evidence.trans.matrix = evidence.trans.matrix[, pos.layer, drop = F],
                                               basic.weights = tdatasets$weights_total [ !is.na(tdatasets$weights_total) & pos.layer ],
                                               main = tlayers$evidence_layer[i],
                                               replace.NA.0 = TRUE, 
                                               p.value.permutation = p.value.permutation, 
                                               n.permutations = n.permutations,
                                               permute.option = permute.option,
                                               # If permute option is "rows per layer", this argument indicates the columns in each layer.
                                               permute.col.groups = tlayers$evidence_layer[i] )
        
        
        
        # Add the probability vector to the table of tci
        ttci.probability = cbind(ttci.probability, 
                                 res.basic.tci.score$gene.tci.probability)
        
        
        
      }
      
      # Compute the P(S>s) to the total tci
      cat(" ++ Computing P(S>s) for total tci. \n")
      
      # Use the whole evidence matrix to permute the scores
      res.basic.tci.score = basic_tci_score( evidence.trans.matrix = evidence.trans.matrix,
                                             basic.weights = tdatasets$weights_total [ !is.na(tdatasets$weights_total) ],
                                             main = "Total TCI",
                                             replace.NA.0 = TRUE, 
                                             p.value.permutation = p.value.permutation, 
                                             n.permutations = n.permutations,
                                             permute.option = permute.option,
                                             # If permute option is "rows per layer", this argument indicates the columns in each layer.
                                             permute.col.groups = tdatasets$evidence_layer [ !is.na(tdatasets$weights_total) ] )
      
      
       
      # Add the probability vector to the table of tci
      ttci.probability = cbind(ttci.probability, 
                                res.basic.tci.score$gene.tci.probability)
      
      # Add colnames to the table and transform to data frame
      colnames(ttci.probability) <- c("ensembl_gene_id",  paste("probability_tci_", tlayers$evidence_layer, sep = ""), "probability_tci_total") # paste("probability_tci_", levidence.layers, sep = "")
      ttci.probability = data.frame(ttci.probability, stringsAsFactors = F)
      print( head(ttci.probability) )
      
    }
    
   
    
    # Return results
    return( list(ttci.score = ttci.score,
                 ttci.probability = ttci.probability,
                 tconfidence.score = tconfidence.score,
                 evidence.trans.matrix = evidence.trans.matrix,
                 ltrans.datasets = ltrans.datasets))
    
    
  }
  
  
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Function to compute the tci for evidence of association using the basic algorithm
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  
  basic_tci_score = function( evidence.trans.matrix = NULL, # Numeric matrix of scaled [0,1] datasets
                              basic.weights = NULL, # Vector of weights for each column of evidence.trans.matrix
                              main = "TCI", # Title for the plots
                              replace.NA.0 = FALSE, # Remplace NAs for 0 in the input matrix
                              p.value.permutation = FALSE, # Compute P(S>s) for score using permutations
                              n.permutations = 10000, # Number of permutations for P(S>s)
                              permute.option = c("both", "columns", "rows", "rows per layer", "rows per column"), # Option to permute rows, columns or both in evidence.trans.matrix
                              permute.col.groups = NULL # If permute option is "rows per layer", this argument indicates the columns in each layer.
    
  ){
    
    
    # Replace NAs for 0's to compute the sum product
    if( replace.NA.0 == TRUE ){
      evidence.trans.matrix [ is.na(evidence.trans.matrix) ] <- 0
    }
    
    # Compute the weighted sum of scores for this evidence matrix
    cat(" ++ Compute the weighted sum of scores for this evidence matrix \n")
    lrow.names.evidence.trans.matrix = rownames(evidence.trans.matrix)
    evidence.trans.matrix = data.matrix(evidence.trans.matrix)
    
    
    print(head(evidence.trans.matrix))
    
    if( ncol(evidence.trans.matrix) > 1 ){
      
      gene.tci.score = round( evidence.trans.matrix %*% as.numeric(basic.weights), 3)  
      
    }else{
      
      gene.tci.score = round( as.numeric(evidence.trans.matrix[,1]) * as.numeric(basic.weights), 3)  
      
    }
    
    
    # Compute P(S>s) of the score using permuations of the data. Please not that:
    # Columns between layers might not permutable because the models generating the data are different
    # Assume the columns are exchangable
    
    # Run the permutation even if only one column in the evidence matrix
    if( p.value.permutation == TRUE & ncol(evidence.trans.matrix) >= 1 ){
      
      
      # Initialise the vector of number of times the score is higher thann observed 
      n.higher.score = rep(0, nrow(evidence.trans.matrix))
      
      # Print type of permutations requested for evidence transformed matrix
      cat(" ++ Running permutations of", permute.option, "\n")
      
      # Loop across permutations
      for( i in 1:n.permutations){
        
        
        # Generate random input matrix by permutating the columns, rows,  both or rows per layer
        if( permute.option == "columns" ){
          
          evidence.trans.matrix.perm = evidence.trans.matrix [, sample(1:ncol(evidence.trans.matrix), size = ncol(evidence.trans.matrix), replace = T)]
          
        }else if( permute.option == "rows" ){
          
          evidence.trans.matrix.perm = evidence.trans.matrix [ sample(1:nrow(evidence.trans.matrix), size = nrow(evidence.trans.matrix), replace = T) ,]
          
        }else if( permute.option == "rows and columns" ){
          
          evidence.trans.matrix.perm = evidence.trans.matrix [ sample(1:nrow(evidence.trans.matrix), size = nrow(evidence.trans.matrix), replace = T), sample(1:ncol(evidence.trans.matrix), size = ncol(evidence.trans.matrix), replace = T)]
        
        }else if( permute.option == "rows per column" ){  
          
          # Initialise permuted matrix
          evidence.trans.matrix.perm <- evidence.trans.matrix 
          
          # Shuffle the rows at each column
          for(j in 1:ncol(evidence.trans.matrix.perm)){
            
            evidence.trans.matrix.perm[,j] <-  evidence.trans.matrix.perm[ sample(1:nrow(evidence.trans.matrix), size = nrow(evidence.trans.matrix), replace = T), j, drop = F]
            
          }
          
        }else if( permute.option == "rows per layer" & !is.null(permute.col.groups) ){
          
          # Identify groups of columns (datasets within each layer) 
          l.groups = unique(permute.col.groups)
          evidence.trans.matrix.perm = NULL
          
          # For each group of columns the rows are permuted in evidence.trans.matrix
          for(j in 1:length(l.groups)){
            
            evidence.trans.matrix.perm = cbind( evidence.trans.matrix.perm,
                                                evidence.trans.matrix [ sample(1:nrow(evidence.trans.matrix), size = nrow(evidence.trans.matrix), replace = T) ,
                                                                        permute.col.groups == l.groups[j] ] 
                                              )
              
          } # End loop across groups of columns
        
        } # End if type of permutation
        
        # Compute the weighted sum of scores for this permuted evidence matrix
        gene.tci.score.perm = data.matrix(evidence.trans.matrix.perm) %*% basic.weights
        
        # Count if the score is higher than observed 
        pos.higher = gene.tci.score.perm > gene.tci.score
        n.higher.score[ pos.higher ] <- n.higher.score[ pos.higher ] + 1
         
        
      } # End loop across permutations
      
   
      # Compute P(S>s)
      gene.tci.probability = n.higher.score / n.permutations
      
      # When no higher score observed, impute based on the number of permutations
      gene.tci.probability[ gene.tci.probability == 0 ] <- 1/n.permutations
      # hist(gene.tci.probability, col = "gray", xlab = "P(S>s)", main = main, ylab = "Number of genes", las = 0, xlim = c(0,1) )
      
      #+ cat( "Summary of P(S>s):")
      #+ print( summary(gene.tci.probability) )
      
    }else{ # If not computing P(S>s) from permutations
      
      gene.tci.probability = NA
      
    }
   
    
    # Return results
    return( list( gene.tci.score = gene.tci.score,
                  gene.tci.probability = gene.tci.probability) )
    
    
    
  }
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Function for finding median differences with bootstrapping
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # Source: Howell. Resampling and Nonparametric Approaches to Data. Chapter 18. pp 12 
  
  bootstrapping.median.differences = function( sample.A = NULL,
                                               sample.B = NULL,
                                               bootstrap.size = 10000,
                                               alternative = c("two.sided", "greater", "less")){
    
    # Find sample size
    sample.A = as.numeric(sample.A)
    sample.B = as.numeric(sample.B)
    n.sample.A = length(sample.A)
    n.sample.B = length(sample.B)
    sample.merged = c(sample.A, sample.B)
    n.sample.merged = length(sample.merged)
    
    # Declare arrays for bootstrapping results (medians and median differences)
    lmedian.A = lmedian.B = lmedian.diff = rep(NA, bootstrap.size)
    
    # Loop for bootstrap
    for( i in 1:bootstrap.size){
      
      # Shuffle samples
      permutated.positions = sample(x = 1:n.sample.merged, 
                                    size = n.sample.merged, 
                                    replace = F)
      
      sample.merged.permutated = sample.merged[ permutated.positions ]
      sample.A.per = sample.merged.permutated[ 1:n.sample.A  ]
      sample.B.per = sample.merged.permutated[ (n.sample.A+1):length(permutated.positions) ]
      
      # Compute median and median differences
      lmedian.A [i] = median(sample.A.per)
      lmedian.B [i] = median(sample.B.per)
      lmedian.diff [i] = lmedian.A [i] - lmedian.B [i] 
      
    }
    
    # Compute observed median difference
    obs.median.diff = median(sample.A) - median(sample.B)
    
    
    # Compute number of permutations were the median was greater or lower than the observed (absolute value)
    # hist(lmedian.diff)
    if( alternative == "greater" ){
      
      count.perm = sum( lmedian.diff > obs.median.diff ) # sum( lmedian.diff > abs(obs.median.diff) ) + sum( lmedian.diff < -abs(obs.median.diff) ) 
      
    }else if( alternative == "less" ){
      
      count.perm = sum( lmedian.diff < obs.median.diff ) # sum( lmedian.diff > abs(obs.median.diff) ) + sum( lmedian.diff < -abs(obs.median.diff) ) 
      
    }else if( alternative == "two.sided"){
      
      count.perm = sum( lmedian.diff > abs(obs.median.diff) ) + sum( lmedian.diff < -abs(obs.median.diff) ) 
      
    }
    
    # Compute p-value # t.test(sample.A, sample.B, alternative = "less" )
    p.value = count.perm / bootstrap.size
    
    # Correct p-value by bootstrap.size
    p.value = max(p.value, 1/bootstrap.size)
    
    return(p.value = p.value)
    
  }
  