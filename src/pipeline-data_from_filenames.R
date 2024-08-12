
### Get files based on data file names:
# (1) .txt file to rename links to data files (replace lane identifiers by 1-3 for consistency across samples)
# NOTE: never rename the raw data files, instead create symbolic links and rename the links
# (2) .csv file for python dictionary
### for


### -------------- load required packages & set path -----------------------------

library(stringr)
p.wd <- "/shared/slate_group1/Shared/methylated_soay/soay_wgbs_pilot_mar2023"
setwd(p.wd)

### -------------- load & format data ----------------------------------------------

FileNames <- list.files(path="Trimmed_renamed", pattern = "*fastq.gz$")

# set up loop to get data from file names:
NewData <- NULL
for(f in seq(1, length(FileNames), by=1)) {
  Parts <- str_split_fixed(FileNames[f],"_",4) # this gets, sample ref numbers (Liverpool & lab), lane identifier (date), lane identifier (lane name), and read pair member tag (or R0) together file file extension
  Data <- data.frame(File=FileNames[f], Parts)
  NewData <- rbind(NewData, Data)
}

# prepare (1) txt file to rename data files (lane)
# get replacement 'lane' identifiers
help <- rep(1:3, each=3)
NewData$lane_dummy <- rep(help, length(unique(NewData$X1))) # file names are sorted by sample, hence repeats of 1,2,3 can be used as order of lanes does not matter 

# create file with old name of links (NewData$File) and new name of links (new.help)
new.help <- paste(NewData$X1, NewData$lane_dummy, NewData$X4, sep="_")
FileRename <- data.frame(original=paste("Trimmed_renamed", NewData$File, sep="/"), new=paste("Trimmed_renamed_2", new.help, sep="/"))

# save file
write.table(FileRename, "src/FileRename.txt", col.names=F, row.names=F, quote=F, sep=" ")


# prepare (2) .csv file for python dictionary
# the file is used to ad read groups to alignment files and requires the lane identifier, flow cell and sample/library name

new.help2 <- paste(NewData$X1, NewData$lane_dummy, sep="_") # new link name / sample name
Flowcell <- NewData$X2
Lane <- str_split_fixed(NewData$X3,"00",2)[,2]
NewData <- data.frame(Sample=new.help2, Flowcell=Flowcell, Lane=Lane)

# reshape data such that dictionary can be created with python
DataOut <- reshape(NewData, idvar="Sample", varying=names(NewData)[-1], v.name="Value", times=names(NewData)[-1], new.row.names=1:(length(names(NewData)[-1])*nrow(NewData)), direction="long")
names(DataOut)[2] <- "Key"
DataOut.temp <- DataOut[!duplicated(DataOut[,1:2]),]

# save file
write.table(DataOut.temp, "src/ReadGroups.help.new.csv", quote=F, sep=",", row.names=F)

