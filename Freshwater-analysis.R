rm(list=ls())
taxa<- read.csv('data/Taxa/INV_OPEN_DATA_TAXA.csv', header=T, sep = ',')
site<- read.csv('data/Taxa/INV_OPEN_DATA_SITE.csv', header=T, sep = ',')
taxoninfo <- read.csv('data/Taxa/OPEN_DATA_TAXON_INFO.csv', header=T, sep = ',')
metrics<- read.csv('data/Taxa/INV_OPEN_DATA_METRICS.csv', header=T, sep=',')
#create summaries for the data to see what your columns and factors are
sumsite<- as.data.frame(summary (site))
summetric<- as.data.frame(summary(metrics))
sumsite

#create a secondary site dataframe to not mess with the original with only columns 
#that have more than 10,000 entries (as the data set is
#30453, set the N/A threshold to 0.3)
combined = site[,!sapply(site, function(x) mean(is.na(x)))>0.7]

#create a copy of the taxa table to manipulate
taxacopy<- taxa
#copy over the column from taxon info, but only the ones that match in both tables in the key column
taxacopy$TAXON_NAME <- taxoninfo$TAXON_NAME[match(taxacopy$TAXON_LIST_ITEM_KEY,
                                                  taxoninfo$TAXON_LIST_ITEM_KEY)]

#make a filtered table with the NA to check for missing taxon names
na<-is.na(taxacopy$TAXON_NAME)
#then use unique to see for anything other than false
unique(na)
remove(na)

#Clean and map the column - UNDER CONSTRUCTION 
metrics <- metrics %>%
  mutate(metrics$REPLICATE_CODE, case_when(
    metrics$REPLICATE_CODE %in% c("Rep1", "Rep 1", "rep 1", "rep1",'Replicate','replicate 1', 'Replicate 1',
                                  'REPLICATE 1', 'rep') ~ "rep1",
    TRUE ~ metrics$REPLICATE_CODE
  )) %>%
  mutate(metrics$REPLICATE_CODE<- as.factor(metrics$REPLICATE_CODE))          # Convert back to factor
unique(metrics$REPLICATE_CODE)

#create the model
mymod<- lmer(COUNT_OF_SAMPLES ~ ALKALINITY +(1|REPORTING_AREA)+(1|BASE_DATA_DATE)+ WIDTH + DEPTH, data=site2)
summary (mymod)
