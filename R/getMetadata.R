## create df for all metrics - all samples generated for Touchstone

basePath = '/mnt/scratch1/Touchstone_data/new_data/'
df_samples <- data.frame(sample_id = dir(basePath),
                         platform = NA, expMat = NA, tx_file = NA, cell_meta = NA)

# add platform based on nomenclature of samples
df_samples$platform[grep('*XR*|*XRS*', df_samples$sample_id)] <- 'Xenium'
df_samples$platform[grep('_CR_', df_samples$sample_id)] <- 'CosMx'
df_samples$platform[c(1:14)] <- c('CosMx', 'CosMx', 'CosMx', 'CosMx', rep('Xenium',10) )
df_samples <- df_samples[-7,] ## corrupted sample
## removing non updated samples from now
df_samples <- df_samples[-c(1:14), ]




# Function to find the first file matching the pattern in each folder and return its path
findFilePath <- function(folderName, pattern) {
  # Construct the full path to the folder
  fullPathToFolder <- paste0(basePath, folderName)

  # Find files in the folder matching the pattern
  files <- list.files(path = fullPathToFolder, pattern = pattern, full.names = TRUE)

  # Return the first matching file's full path, or NA if no file was found
  if (length(files) > 0) {
    return(files[1])
  } else {
    return(NA)
  }
}

# Apply the function to each folder name and store the results in a new column
df_samples$expMat <- sapply(df_samples$sample_id, function(x) { findFilePath(folderName = x, pattern = "*exprMat*|cell_feature_matrix") } )
df_samples$tx_file <- sapply(df_samples$sample_id, function(x) { findFilePath(folderName = x, pattern = "*tx_unique*|*tx_file_unique*|transcripts") } )
df_samples$cell_meta <- sapply(df_samples$sample_id, function(x) { findFilePath(folderName = x, pattern = "*cells.csv.gz*|metadata") } )
head(df_samples)
