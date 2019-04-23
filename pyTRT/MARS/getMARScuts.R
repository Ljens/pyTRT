# Script is a standalone MARS function to be called from python
# Returns the hinge locations, or cut points


# Get args from command line
# data_file should be a csv with two columns:
  # Time in hours, column called "Hours"
  # Average fluid temp, column call "Temp_ave"

# Input test:
#data_file = "/Users/Linden/Google Drive/PHD/Models/Jnotebooks/MARS/to_earth.csv"


#temp = measured temperature
#nk = max terms, general 10

args = commandArgs(trailingOnly=TRUE)
data_file = args[1]
nk = args[2]
#cat('###############')
#cat(nk)
#cat('###############')
  
# Use read.csv to turn the datafile input into a dataframe
# NOTE: if data_file is an actual file, use read.csv(data_file, ...)
#       if data_file is a text string in csv format straight from python, use read.csv(text = data_file, ...)
#data = read.csv(text = data_file,  header=TRUE, stringsAsFactors=FALSE, sep=",", colClasses=("numeric"))
data = read.csv(data_file,  header=TRUE, stringsAsFactors=FALSE, sep=",", colClasses=("numeric"))

data$logHours = log(data$Hours)

# Import the earth package
library(mda)

mars_fit=mars(data$logHours,data$Temp_ave,degree=1,nk=as.numeric(nk)) #Fit a MARS model to the data
cuts=c(min(data$logHours),sort(unique(mars_fit$cuts)[-1])) #Take the unique cut values

# Print the cut points in log space to the commandline
cat(cuts)
#cat(mars_fit$fitted)

