library(Revelio) # For Revelio (HeLa cells) data and method. 
library(dplyr) # For groupby
library(gmodels) # For fast.prcomp (works well when num_rows << num_cols)
library(data.table) # For fwrite(matrix_name, file_name)
# setwd("~/Desktop/Research/structured-priors/experiments/cell-cycle/")

# Store in a Revelio object the HeLa Cell dataset used in Appendix of 
# the original Revelio paper and run PCA as well as Revelio on it, 
# saving the results in the object.
extract_data_results <- function() {
    my_data <- createRevelioObject(rawData = revelioTestData_rawDataMatrix,
                                  cyclicGenes = revelioTestData_cyclicGenes)
    my_data <- getCellCyclePhaseAssignInformation(dataList = my_data) # need to run this to obtain their inferred cell phase assignments
    my_data <- getPCAData(dataList = my_data)
    my_data <- getOptimalRotation(dataList = my_data)
}

my_data <- extract_data_results()

pc_data <- my_data@transformedData$pca$data 
fwrite(pc_data, "../data/pca-hela-cells.csv")

count_data <- myData@DGEs$countData[myData@geneInfo$variableGenes,] # subset it for variable genes 
fwrite(count_data, "../data/hela-X-variable-genes.csv")
count_data <- cbind(data.frame(count_data), cell_phase_assignments) # add cell phase assignements for easy groupby
colnames(count_data)[1346] <- "cell_phases"

scaled_count_data <- myData@DGEs$scaledData[myData@geneInfo$variableGenes,]
fwrite(scaled_count_data, "../data/hela-X-variable-genes-scaled.csv", row.names=TRUE)
fwrite(myData@DGEs$scaledData, "../data/hela-X-all-genes-scaled.csv", row.names=TRUE)

cell_phase_assignments <- myData@cellInfo$ccPhase
fwrite(as.vector(cell_phase_assignments), "hela-cell-phase-assignments.csv")