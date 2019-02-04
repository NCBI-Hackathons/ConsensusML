#Jenny Smith

#Purpose: Concatenate the expression data-sets downloaded from TCGA/TARGET from GDC or any patient level data
#eg. each individual patient has a single expression-file 


catExpnData <- function(filenames,regex, cols, header=FALSE,removeFirstLine=FALSE, sep="\t"){
  library(magrittr)
  options(stringsAsFactors = FALSE)
  #filenames is a character vector of all filenames. 
  #regex is a string with the pattern to extract the patient ID , eg "^.+(Kasumi|MV4)", from filenames 
  #cols is the character vector or numeric vector of the columns to select and concatenate. 
  
  extract_cols <-function(filename,cols,rmFirstLine=FALSE){
    
    if(all(rmFirstLine & header)){
      aFile <- readLines(filename)[-1] #remove first line with extra info. 
      aFile <- str_split_fixed(aFile, pattern = "\t",n = length(cols)) %>% #split into a matrix
        set_colnames(.[1,] ) %>%  #set colnames from the first line 
        .[-1, ] #remove the header row from matrix
    }else{
      aFile <- read.delim(filename, sep=sep, header=header, as.is=TRUE)
    }
    
    output <- list()
    for ( k in 1:length(cols)){
      colname <- cols[k]
      col <- aFile[,colname]
      output[[colname]] <- col
    }
    return(output)
  }
  
  combineColumns <- function(extract_cols.res,colname){
    sapply(extract_cols.res, '[[', colname)
  }
  
  
  IDs <- gsub(regex, "\\1", filenames)
  
  columns <- lapply(filenames,extract_cols,cols=cols, rmFirstLine=removeFirstLine) %>%
    set_names(IDs)
  
  catedMatrices <- lapply(cols, combineColumns, extract_cols.res=columns)  %>%
    set_names(cols)
  
  
  return(catedMatrices)
}



