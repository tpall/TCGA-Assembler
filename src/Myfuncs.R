
# convert .RData -> .rdb/.rdx
tools:::makeLazyLoadDB(local({load(travRes); environment()}), 
                       "DirectoryTraverseResult_5-Jun-2015")
lazyLoad("DirectoryTraverseResult_5-Jun-2015")
ls()
x 

MyDownloadMethylationData <- function(traverseResultFile, saveFolderName, cancerType, assayPlatform, tissueType = NULL, inputPatientIDs = NULL, outputFileName = "")
  library(stringr)
  {
  options(warn=-1)
  
  writeLines("**********************************************************************************")
  writeLines("");
  writeLines(paste("Download DNA methylation data of ", cancerType, " patients.", sep = ""))
  
  # Check whether specific TCGA patient IDs are inputted. 
  if (!is.null(inputPatientIDs))
  {
    inputPatientIDs = str_to_upper(gsub(pattern = "\\.", replacement = "-", x = inputPatientIDs),"en")
  }
  
  if ((outputFileName != "") & (!is.null(outputFileName)))
  {
    outputFileName = paste(outputFileName, "__", sep = "")
  }
  
  # load directory traverse result and create folder to store output files
  writeLines("Load information of TCGA data files.")
  load(traverseResultFile)
  dir.create(path = saveFolderName, recursive = TRUE)
  SpecificID = grep(pattern = str_to_upper(paste("/", cancerType, "/cgcc/jhu-usc\\.edu/", assayPlatform, "/", sep = ""),"en"), x = upper_file_url, ignore.case = FALSE);
  MeLevel3ID = SpecificID[grep(pattern = str_to_upper("Level_3","en"), x = upper_file_url[SpecificID], ignore.case = FALSE)];
  
  # search for the Sample and Data Relationship Format (SDRF) file of the specified platform and cancer type
  ind = SpecificID[grepEnd(pattern = str_to_upper("\\.sdrf\\.txt","en"), x = upper_file_url[SpecificID], ignore.case = FALSE)];
  if (length(ind) == 0)
  {
    writeLines("Program exited due to missing SDRF file."); 
    return();
  }
  if (length(ind) > 1)
  {
    URL = GetNewestURL(AllURL = file_url[ind]);
  }else{
    URL = file_url[ind];
  }
  downloadResult = urlReadTable(url = URL);
  if (downloadResult$errorFlag != 0)
  {
    writeLines("Error in downloading SDRF file.");
    return();
  }
  sdrf = str_to_upper(downloadResult$data,"en");
  
  # Process SDRF file, identify the columns of level 3 data file name and TCGA sample barcodes.
  level_3_filename_column = max(grep(pattern = "Data Matrix File", x = sdrf[1, ], ignore.case = TRUE));
  DataLevelColID = max(grep(pattern = "TCGA Data Level", x = sdrf[1, ], ignore.case = TRUE));
  TCGABarcodeID = min(grep(pattern = "TCGA Barcode", x = sdrf[1, ], ignore.case = TRUE));
  colnames(sdrf) = sdrf[1, ];
  sdrf = unique(sdrf[2:dim(sdrf)[1], , drop = FALSE]);
  sdrf = sdrf[!duplicated(sdrf[, level_3_filename_column]), c(TCGABarcodeID, DataLevelColID, level_3_filename_column), drop = FALSE];
  SDRFID = sort(union(which(sdrf[, 2] == "LEVEL_3"), which(sdrf[, 2] == "LEVEL 3")), decreasing = FALSE);
  if (length(SDRFID) == 0)
  {
    writeLines("Error: there are no Level 3 data");
    return();
  }
  sdrf = sdrf[SDRFID, , drop = FALSE];  
  
  # If specific patient TCGA barcodes are inputted, only download the specified samples.
  if (!is.null(inputPatientIDs))
  {
    indInputPatientID = c();
    for (i in 1:length(inputPatientIDs))
    {
      indInputPatientID = c(indInputPatientID, grepBeginning(pattern = inputPatientIDs[i], x = sdrf[, 1], ignore.case = FALSE));
    }
    if (length(indInputPatientID) == 0)
    {
      writeLines("No Level 3 data for the inputted TCGA barcodes.");
      return();      
    }else{
      sdrf = sdrf[indInputPatientID, , drop = FALSE];
    }
  }
  
  # Download data of specified tissue
  if (!is.null(tissueType))
  {
    SampleType = cbind(Options = c("TP", "TR", "TB", "TRBM", "TAP", "TM", "TAM", "THOC", "TBM", "NB", "NT", "NBC", "NEBV", "NBM"),
                       Code = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14"));  
    sdrf = sdrf[substr(sdrf[, 1], 14, 15) %in% SampleType[SampleType[, "Options"] %in% tissueType, "Code"], , drop = FALSE];
  }
  
  if (dim(sdrf)[1] == 0)
  {
    writeLines("No available data.");
    return();
  }
  
  # Download data files of all samples.
  left_columns = NULL;
  AllPosition = NULL;
  exp_names = NULL;
  data = NULL;
  for (i in 1:dim(sdrf)[1])
  {
    time1 = proc.time();
    sample_TCGA_id = sdrf[i, 1];
    ind = MeLevel3ID[grepEnd(pattern = str_to_upper(sdrf[i, 3],"en"), x = upper_file_url[MeLevel3ID], ignore.case = FALSE)];
    if (length(ind) == 0)
    {
      next;
    } 
    if (length(ind) > 1)
    {
      URL = GetNewestURL(AllURL = file_url[ind]);
    }else{
      URL = file_url[ind];
    }
    downloadResult = urlReadTable(url = URL);
    if (downloadResult$errorFlag != 0)
    {
      next;
    }
    s = downloadResult$data;
    s = s[2:dim(s)[1], , drop = FALSE];
    
    chr = rep(0, dim(s)[1]);
    for (j in 1:22)
    {
      IDj = which(s[, 4] == as.character(j));
      chr[IDj] = j;
    }
    IDj = which(str_to_upper(s[, 4],"en") == "X");
    chr[IDj] = 23;      
    IDj = which(str_to_upper(s[, 4],"en") == "Y");
    chr[IDj] = 24; 
    position = rep(0, dim(s)[1]);
    position[2:length(position)] = as.numeric(s[2:dim(s)[1], 5]);
    Yj = chr*(10e+10) + position;
    orderIDj = order(Yj, decreasing = FALSE);
    Yj = Yj[orderIDj];
    s = s[orderIDj, , drop = FALSE];
    
    # Need to check whether every data file has the same methylation probes
    if (is.null(left_columns))
    {
      left_columns = s[, c(1, 3, 4, 5), drop = FALSE];
      AllPosition = Yj;
      exp_names = sample_TCGA_id;
      data = s[2:dim(s)[1], 2, drop = FALSE];
    }else{
      if (sum(AllPosition != Yj) > 0)
      {
        next;
      }
      exp_names = c(exp_names, sample_TCGA_id);
      data = cbind(data, s[2:dim(s)[1], 2, drop = FALSE]);
    }
    
    time = proc.time() - time1;
    writeLines(paste("Downloaded - ", cancerType, "__jhu-usc.edu__", assayPlatform, " - sample ", i, " out of ", dim(sdrf)[1], ". ", round(time[3], digits = 1), " seconds elapsed.", sep = ""));
  }
  
  writeLines("Save data to local disk.");
  ID = str_locate_all(traverseResultFile, "_")[[1]];
  ID = ID[dim(ID)[1], 2];
  filename = paste(saveFolderName, "/", outputFileName, cancerType, "__jhu-usc.edu__", assayPlatform, "__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");  
  header = c(left_columns[1, ], exp_names);
  left_columns = left_columns[2:dim(left_columns)[1], , drop = FALSE];
  ID = grepBeginning(pattern = "NA", x = left_columns[, 3], ignore.case = TRUE)  
  ID = sort(setdiff(1:dim(left_columns)[1], ID), decreasing = FALSE);
  left_columns = left_columns[ID, , drop = FALSE];
  data = data[ID, , drop = FALSE];
  write.table(rbind(header, cbind(left_columns, data)), file = filename, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE);  
  
  writeLines("");
  writeLines("**********************************************************************************");
  writeLines("\n");
  
  options(warn=0);
  
  # Return downloaded data
  downloadedData = rbind(header, cbind(left_columns, data));
  rownames(downloadedData) = NULL;
  colnames(downloadedData) = NULL;
  downloadedData;
}