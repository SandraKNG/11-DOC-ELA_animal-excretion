#setwd("D:/Xenopoulos Lab/PARAFAC/PARAFAC 2021 documents/Test for lab")

# 1. Make a list of the .csv filenames:
file_list = list.files(pattern="*.csv")
write.csv(file_list,"FinalChangeLog.csv")

# 2. Edit FinalChangeLog.csv to have the following columns (start by deleting first column):
## a) OLDNAME (whatever.csv)
## b) NEWNAME (PARAFAC_001, PARAFAC_002, ..., PARAFAC_00n)

# 3. Read in edited changelog:

y<-read.csv("FinalChangeLog.csv")

# 4. Rename EEM files:
data = NULL
for (j in 1:nrow(y)){
  data<-read.csv(paste0(y[j,1],""),header=FALSE)
  write.table(data,paste0(y[j,2],".csv"),row.names=FALSE,col.names=FALSE, sep = ",")
}

# 5. Transfer the new PARAFAC_### EEMs to their own folder CONTAINING NO OTHER .CSV FILES! 
#    This is the folder you will use in MATLAB to make your big matrix. If there are any
#    other .csv files in the folder they will mess things up.