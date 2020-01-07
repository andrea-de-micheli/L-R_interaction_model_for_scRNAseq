#Ligand-receptor interaction model for single-cell RNA-sequencing experiments.
#The code is provided as an example without any warranty.
#Copyright © 2020 Andrea De Micheli

###############################
#                             #
#    LIGAND-RECEPTOR (L-R)    #
#    database initialisation  #
#                             #  
###############################

#STEP 1 -- Import ligand-receptor database and gene expression dataset (gene expression matrix or Seurat object):
receptor_ligand_dataset = read.table("./database.txt", header = T, sep = "\t")
rownames(receptor_ligand_dataset) = receptor_ligand_dataset$pair #Make sure the row names are the paris and that ligands and receptor columns are respectively called ligand_symbol and receptor_symbol
expression_dataset = load("expression_matrix.RData")

#STEP 2 -- Subset database from genes in the experiment expression dataset:
#This function will subset the L-R database from the genes in the gene expression matrix (experiment).
#Will return a list where [[1]] are the list of pairs in the experiment, and [[2]] the subsetted RL database
RLSubset = function(receptor_ligand_dataset, dataset) {
  
  list_RLpair_remove = c()
  for (i in 1:length(receptor_ligand_dataset$pair)) {
    
    if (identical(intersect(as.character(receptor_ligand_dataset$ligand_symbol[i]),rownames(dataset)), character(0)) || 
        identical(intersect(as.character(receptor_ligand_dataset$receptor_symbol[i]),rownames(dataset)), character(0))) {
      
      list_RLpair_remove = c(list_RLpair_remove, as.character(receptor_ligand_dataset$pair[i]))
    }
    
  }
  
  receptor_ligand_dataset_subset = receptor_ligand_dataset[!(rownames(receptor_ligand_dataset) %in% list_RLpair_remove),]
  list_RLpair = rownames(receptor_ligand_dataset_subset)
  
  return(list(list_RLpair,receptor_ligand_dataset_subset))
  
}

#Subsetting the receptor-ligand dataset with only the receptors and ligands expressed in experiment:
receptor_ligand_dataset = RLSubset(receptor_ligand_dataset, expression_dataset)
list_RLpair = receptor_ligand_dataset[[1]]

#Extract the list of expressed ligand and receptor genes:
list_Receptor = unique(as.vector(receptor_ligand_dataset[[2]]$receptor_symbol))
list_Ligand = unique(as.vector(receptor_ligand_dataset[[2]]$ligand_symbol))
list_RLgenes = unique(c(list_Ligand,list_Receptor))

################################
#                              #
#   Gene expression averages   #
#     and differentially       #
#     expressed receptors      #
#                              #  
################################

#STEP 3 -- Calculate gene expression average across labeled cell populations
#The example probided below is for a Seurat v3 (Satija Lab) dataset:

#Calculating average expression values using raw counts:
require(Seurat)
expression_dataset = SetIdent(expression_dataset, value = "cell_annotation")
expression_dataset_averages = AverageExpression(expression_dataset, genes.use = list_RLgenes, slot = "counts") #Non-normalized expression values (counts)

#STEP 4 -- Obtain list of differentially expressed receptors (or ligands)
DE_receptors = read.table("./DE_gene_list.txt")

################################
#                              #
#     Compute L-R matrix       #
#                              #  
################################

#The algorithm  will take the average expression matrix, the subsetted L-R database, the list of differentially expressed receptor or ligand, and the target cell type to either:
#1. Calculate interactions between DE receptors on "type" cell type vs. ligands expressed on all other cell types (default, contrast = "receptor")
#2. Calculate interactions between DE ligands on "type" cell type vs. receptors expressed on all other cell types (contrast = "ligand")
#See STEP 5 (line 148) for execution in context.
ComputeRLmatrixType = function(matrix, receptor_ligand_dataset, type, diff_genes, contrast = "receptor") {
  
  require(dplyr)
  require(data.table)
  
  RLPairMatrix = data.frame(matrix(ncol = dim(matrix)[2], nrow = dim(receptor_ligand_dataset)[1]))
  colnames(RLPairMatrix) = colnames(matrix)
  rownames(RLPairMatrix) = rownames(receptor_ligand_dataset)
  
  #If contrasting receptors on one cell type to ligands expressed by other cell types:
  if (contrast == "receptor") {
    
    for (pair in 1:dim(receptor_ligand_dataset)[1]) {
      
      #Receptor expression value for target cell type
      receptor = matrix[as.character(receptor_ligand_dataset$receptor_symbol[pair]),type]
      
      #If receptor value is not zero and if the gene is differentially expressed  
      if (receptor != 0 && 
          any(diff_genes %in% receptor_ligand_dataset$receptor_symbol[pair]) == TRUE) {  
        
        #For each other cell type
        for (i in 1:dim(matrix)[2]) {
          
          ligand = matrix[as.character(receptor_ligand_dataset$ligand_symbol[pair]),i]
          RLPairMatrix[rownames(receptor_ligand_dataset)[pair], colnames(matrix)[i]] = ligand * receptor
          
        }
      } 
      else {
        RLPairMatrix[rownames(receptor_ligand_dataset)[pair], ] = NA
      }
    }
  }
  
  #If contrasting ligands on one cell type to receptors expressed by other cell types:
  else if (contrast == "ligand") {
    
    for (pair in 1:dim(receptor_ligand_dataset)[1]) {
      
      #Ligand expression value for target cell type
      ligand = matrix[as.character(receptor_ligand_dataset$ligand_symbol[pair]),type]
      
      #If ligand value is not zero and if the gene is differentially expressed  
      if (ligand != 0 && 
          any(diff_genes %in% receptor_ligand_dataset$ligand_symbol[pair]) == TRUE) {  
        
        #For each other cell type
        for (i in 1:dim(matrix)[2]) {
          
          receptor = matrix[as.character(receptor_ligand_dataset$receptor_symbol[pair]),i]
          RLPairMatrix[rownames(receptor_ligand_dataset)[pair], colnames(matrix)[i]] = receptor * ligand
          
        }
      } 
      else {
        RLPairMatrix[rownames(receptor_ligand_dataset)[pair], ] = NA
      }
    }
    
  }
  
  #Removing NAs and rows with only 0
  RLPairMatrix = RLPairMatrix[complete.cases(RLPairMatrix),]
  RLPairMatrix = RLPairMatrix[rowSums(RLPairMatrix) > 0,]
  
  return(RLPairMatrix)
  
}

#STEP 5 -- Computation of L-R score matrix:
RLScoreMatrix = ComputeRLmatrixType(expression_dataset_averages, receptor_ligand_dataset[[2]], diff_genes = DE_receptors, type = "the type of cell to consider for receptors", contrast = "receptor")

#STEP 6 -- Example how to plot the score matrix as a heatmap:
require(pheatmap)
require(RColorBrewer)
pheatmap(RLScoreMatrix, scale = "row", cluster_rows = T, cluster_cols = F, main = "title", color = colorRampPalette(colors = c("red", "#FEFEFE", "blue"))(100), border_color = NA, cellwidth = 10, cellheight = 10)
