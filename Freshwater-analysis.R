#Clear the environment and plots
rm(list=ls())
graphics.off()

#Load in or install the necessary packages 
packages <- c("lme4", "tidyverse", "lubridate", "foreign", "data.table",
              "hutilscpp", "sf", "RANN", "geosphere", 'DHARMa', 'ggeffects','patchwork')

# Install any that aren't already installed
missing <- packages[!packages %in% installed.packages()[, "Package"]]
if (length(missing)) install.packages(missing)

# Load them all
invisible(lapply(packages, library, character.only = TRUE))
remove('missing', 'packages')

#Load in the data
taxa<- read.csv('data/Taxa/INV_OPEN_DATA_TAXA.csv', header=T, sep = ',')
site<- read.csv('data/Taxa/INV_OPEN_DATA_SITE.csv', header=T, sep = ',')
taxoninfo <- read.csv('data/Taxa/OPEN_DATA_TAXON_INFO.csv', header=T, sep = ',')
metrics<- read.csv('data/Taxa/INV_OPEN_DATA_METRICS.csv', header=T, sep=',')

#Create a secondary site dataframe to not mess with the original with only columns 
#That have more than 10,000 entries (as the data set is
#30453, set the N/A threshold to 0.3)
combined = site[,!sapply(site, function(x) mean(is.na(x)))>0.7]

#create a copy of the taxa table to manipulate
taxacopy<- taxa
#copy over the column from taxon info, but only the ones that match in both tables in the key column
taxacopy$TAXON_NAME <- taxoninfo$TAXON_NAME[match(taxacopy$TAXON_LIST_ITEM_KEY,
                                                  taxoninfo$TAXON_LIST_ITEM_KEY)]

#Add all of the columns necessary to the table to prepare for the temperature data
taxacopy$DISCHARGE<- combined$DISCHARGE[match(taxacopy$SITE_ID,combined$SITE_ID)]
taxacopy$ALKALINITY<- combined$ALKALINITY[match(taxacopy$SITE_ID,combined$SITE_ID)]
taxacopy$MIN_SAMPLE_DATE<- combined$MIN_SAMPLE_DATE[match(taxacopy$SITE_ID,combined$SITE_ID)]
taxacopy$MAX_SAMPLE_DATE<- combined$MAX_SAMPLE_DATE[match(taxacopy$SITE_ID,combined$SITE_ID)]
taxacopy$MAX_SAMPLE_DATE<- combined$MAX_SAMPLE_DATE[match(taxacopy$SITE_ID,combined$SITE_ID)]
taxacopy$WATER_BODY<- combined$WATER_BODY[match(taxacopy$SITE_ID,combined$SITE_ID)]
taxacopy$BMWP_N_TAXA<- metrics$BMWP_N_TAXA[match(taxacopy$SAMPLE_ID,metrics$SAMPLE_ID)]
taxacopy$BMWP_ASPT<- metrics$BMWP_ASPT[match(taxacopy$SAMPLE_ID,metrics$SAMPLE_ID)]
taxacopy$WHPT_N_TAXA<- metrics$WHPT_N_TAXA[match(taxacopy$SAMPLE_ID,metrics$SAMPLE_ID)]
taxacopy$WHPT_ASPT<- metrics$WHPT_ASPT[match(taxacopy$SAMPLE_ID,metrics$SAMPLE_ID)]
taxacopy$EASTING<- combined$EASTING[match(taxacopy$SITE_ID,combined$SITE_ID)]
taxacopy$NORTHING<- combined$NORTHING[match(taxacopy$SITE_ID,combined$SITE_ID)]
taxacopy$NGR_PREFIX<- combined$NGR_PREFIX[match(taxacopy$SITE_ID,combined$SITE_ID)]

#Filter the table by the date range which matches that of the temperature data in both
#min and max sample date columns
combinedriver<- taxacopy
combinedriver$MIN_SAMPLE_DATE <- dmy(combinedriver$MIN_SAMPLE_DATE)
combinedriver$MIN_SAMPLE_DATE <- ymd(combinedriver$MIN_SAMPLE_DATE)

combinedriver$MAX_SAMPLE_DATE <- dmy(combinedriver$MAX_SAMPLE_DATE)
combinedriver$MAX_SAMPLE_DATE <- ymd(combinedriver$MAX_SAMPLE_DATE)

combinedriver<- filter(combinedriver, combinedriver$MIN_SAMPLE_DATE >= '1982-01-01' &
                                     combinedriver$MIN_SAMPLE_DATE <= '2011-12-31')

combinedriver<- filter(combinedriver, combinedriver$MAX_SAMPLE_DATE >= '1982-01-01' &
                         combinedriver$MAX_SAMPLE_DATE <= '2011-12-31')
#Load in temperature data

segments<- read.dbf('data/WaterTemperature/SEGMENTS.DBF', as.is = F)
tempfiles <- list.files(path = 'data/WaterTemperature/segments/', pattern = '\\.csv$',
                        full.names = TRUE)
#Create a function to put each table as one column with the name of the file as the col name
readfun <- function(file) {
  #Use the file name without the extension as column name
  col_name <- tools::file_path_sans_ext(basename(file))
  #Read the CSV and select the first column
  data <- read_csv(file, col_types = cols(.default = col_double()))
  #Return a tibble with one column named after the file
  tibble(!!col_name := data[[1]])
}

#Read each CSV as a column and combine into one data frame
tempdata <- map_dfc(tempfiles, readfun)
#Remove the first letter so that each column can be matched with the ID in 'segments'
tempdata<- rename_with(tempdata,~ substring(., 2))

#Convert the easting and northing to longitude and latitude
#df2<- combinedriver 
#df2_coords<- df2 %>%
 # st_as_sf(coords = c("EASTING", "NORTHING"), crs = 27700) %>%
  #st_transform(4326) %>%
  #st_coordinates() %>%
  #as.data.frame()

#combinedriver$LAT = df2_coords$Y
#combinedriver$LON = df2_coords$X

# Match the coordinates using the hutilscpp package
#temp_coords<- as.data.frame(subset(segments, select = c(FNODE_LAT, FNODE_LON)))%>%
#st_as_sf(coords = c("FNODE_LON", "FNODE_LAT"), crs = 4326)

#df2_coords <- st_as_sf(df2_coord)%>%
#st_transform(4326)



##Claude method
segments$FNODE_LON<- round(segments$FNODE_LON, 6)
segments$FNODE_LAT<- round(segments$FNODE_LAT, 6)

#Function creation: letter -> 0-based index (A=0)
utf8ToInt2 <- function(ch) utf8ToInt(ch) - utf8ToInt("A")

#Convert a 2-letter NGR prefix to its full BNG easting/northing offset

#Vectorised NGR prefix -> full BNG offset
ngr_offset <- function(prefix) {
  prefix <- toupper(trimws(prefix))
  # 0-based letter index via alphabet lookup (vectorised)
  l1 <- match(substr(prefix, 1, 1), LETTERS) - 1
  l2 <- match(substr(prefix, 2, 2), LETTERS) - 1
  # 'I' is skipped in the grid: shift letters after H (index 7) down by one
  l1 <- ifelse(l1 > 7, l1 - 1, l1)
  l2 <- ifelse(l2 > 7, l2 - 1, l2)
  e100 <- ((l1 - 2) %% 5) * 5 + (l2 %% 5)
  n100 <- (19 - (l1 %/% 5) * 5) - (l2 %/% 5)
  data.frame(e_off = e100 * 100000, n_off = n100 * 100000)
}

#Sanity check
ngr_offset(c("SV", "SX", "TQ", "SO"))

#finish off the easting and northing prefixes
off <- ngr_offset(combinedriver$NGR_PREFIX)

combinedriver$FULL_EASTING  <- off$e_off + combinedriver$EASTING
combinedriver$FULL_NORTHING <- off$n_off + combinedriver$NORTHING

#Transform the full BNG to lat/lon
latlon_to_xyz <- function(lat, lon) {
  lat_r <- lat * pi / 180
  lon_r <- lon * pi / 180
  cbind(
    x = cos(lat_r) * cos(lon_r),
    y = cos(lat_r) * sin(lon_r),
    z = sin(lat_r)
  )
}
pts <- st_as_sf(
  combinedriver,
  coords = c("FULL_EASTING", "FULL_NORTHING"),
  crs = 27700,
  remove = FALSE
)
pts_wgs <- st_transform(pts, crs = 4326)

coords <- st_coordinates(pts_wgs)
combinedriver <- st_drop_geometry(pts_wgs)
combinedriver$LON <- coords[, "X"]
combinedriver$LAT <- coords[, "Y"]

range(combinedriver$FULL_EASTING)   # expect ~80,000–660,000
range(combinedriver$FULL_NORTHING)  # expect ~0–660,000
range(combinedriver$LAT)            # expect ~50–56
range(combinedriver$LON)            # expect ~ -6 to +2

combinedriver$LON<- round(combinedriver$LON, 6)
combinedriver$LAT<- round(combinedriver$LAT, 6)
unique_coords <- combinedriver %>% distinct(LAT, LON)

drv_xyz <- latlon_to_xyz(unique_coords$LAT, unique_coords$LON)
seg_xyz <- latlon_to_xyz(segments$FNODE_LAT, segments$FNODE_LON)

nn <- nn2(data = seg_xyz, query = drv_xyz, k = 1)
matched_idx <- nn$nn.idx[, 1]

lookup <- unique_coords %>%
  bind_cols(segments[matched_idx, ] %>% rename_with(~ paste0(.x, "_seg")))

lookup$distance_m <- distHaversine(
  cbind(lookup$LON, lookup$LAT),
  cbind(lookup$FNODE_LON_seg, lookup$FNODE_LAT_seg)
)

result <- combinedriver %>% left_join(lookup, by = c("LAT", "LON"))

# The moment of truth:
summary(result$distance_m)
sum(is.na(result$distance_m))   # should be 0

#sanity checks
result %>%
  arrange(desc(distance_m)) %>%
  select(LAT, LON, NGR_PREFIX, FNODE_LAT_seg, FNODE_LON_seg, distance_m) %>%
  head(20)

set.seed(1)
plot_sample <- result[sample(nrow(result), 500), ]


##Match ID CREATION

result <- result %>%
  mutate(match_id = as.integer(factor(
    paste(LAT, LON, FNODE_LAT_seg, FNODE_LON_seg)
  )))

#Temperature + Date

#Generate the calendar days from 10/01/1982

tempdata <- tempdata %>%
  mutate(date = seq(as.Date("1982-01-10"),
                    by = "day",
                    length.out = nrow(tempdata)))

##Reshape tempdata to long

##Check necessary segments (i.e. the unique ones) should be about 4840
needed_ids <- unique(as.character(result$OBJECTID_seg))
length(needed_ids)

##Keep only 'date' + the needed segment columns
keep_cols <- c("date", intersect(needed_ids, names(tempdata)))
temp_dt   <- as.data.table(tempdata)[, ..keep_cols]

##Pivot the table from wide to long (rows become columns) by 'melting' the much smaller table
temp_long <- melt(temp_dt,
                  id.vars = "date",
                  variable.name = "object_id",
                  value.name = "temperature")
temp_long[, object_id := as.character(object_id)]
setkey(temp_long, object_id, date)

##Build res_dt and add object_id to store 'result' as a data table
res_dt <- as.data.table(result)
res_dt[, object_id := as.character(OBJECTID_seg)]

#Get the unique combos of segment id, and the start/end of sampling
combos <- unique(res_dt[, .(object_id, MIN_SAMPLE_DATE, MAX_SAMPLE_DATE)])

#Joining the temperature to a new table called avg temp via the unique sets
#Use a non-equi join method to get the average temperature over each date range
avg_temp <- temp_long[combos,
                      on = .(object_id,
                             date >= MIN_SAMPLE_DATE,
                             date <= MAX_SAMPLE_DATE),
                      .(mean_temp = mean(temperature, na.rm = TRUE),
                        n_days    = .N),
                      by = .EACHI]

#Rename every column by its position in one go so they match the rest of the data 
setnames(avg_temp,
         c("object_id", "MIN_SAMPLE_DATE", "MAX_SAMPLE_DATE", "mean_temp", "n_days"))

#Join the tables back onto the full result table
result <- merge(res_dt, avg_temp,
                by = c("object_id", "MIN_SAMPLE_DATE", "MAX_SAMPLE_DATE"),
                all.x = TRUE)

#Final touches
#Ensure it is in Date format
metrics$SAMPLE_DATE <- as.Date(metrics$SAMPLE_DATE, format = "%d/%m/%Y")

#Give one row per SAMPLE_ID + date, then restrict to the temperature series window
sample_dates <- metrics %>%
  distinct(SAMPLE_ID, SAMPLE_DATE) %>%
  filter(SAMPLE_DATE >= as.Date("1982-01-10") &
           SAMPLE_DATE <= as.Date("2011-12-31"))

#Sanity check, Should be 0 (making sure no ID appears more than once)
sample_dates %>% count(SAMPLE_ID) %>% filter(n > 1)

#For extra information: unique samples either side of the date boundary 
metrics %>%
  distinct(SAMPLE_ID, SAMPLE_DATE) %>%
  filter(SAMPLE_DATE < as.Date("1982-01-10") |
           SAMPLE_DATE > as.Date("2011-12-31")) %>%
  nrow()

#Put it back into result
result <- merge(result,
                as.data.table(sample_dates),
                by = "SAMPLE_ID",
                all.x = TRUE)

#Merge temperature again
result <- merge(
  result,
  temp_long,
  by.x = c("object_id", "SAMPLE_DATE"),
  by.y = c("object_id", "date"),
  all.x = TRUE
)
setnames(result, "temperature", "temp_on_day")

#Sanity checks:
sum(is.na(result$temp_on_day))     #Ideally only the rows that had NA SAMPLE_DATE
summary(result$temp_on_day)        #Plausible seasonal river temps
head(result[, .(SAMPLE_ID, object_id, SAMPLE_DATE, temp_on_day, mean_temp)])

#These two should match if all in-range dates found a temperature
sum(is.na(result$temp_on_day))
sum(is.na(result$SAMPLE_DATE))

#Sort to be statistically valid as the structure of result is entirely duplications
#which would pseudoreplicate the response variable
aspt_data <- result %>% distinct(SAMPLE_ID, .keep_all = TRUE)
aspt_model_data <- aspt_data %>%
  filter(!is.na(DISCHARGE)) %>%
  mutate(
    DISCHARGE = as.factor(DISCHARGE),
    temp_z    = as.numeric(scale(temp_on_day)),
    alk_z     = as.numeric(scale(ALKALINITY))
  )

#Filter out the EPT taxa from the taxoninfo and add them to result to prepare for 
#taxa based modelling
ept_taxa <- taxoninfo %>%
  filter(grepl("Ephemeroptera|Plecoptera|Trichoptera", TAXON_GROUP_NAME))
ept_taxa <- ept_taxa %>%
  mutate(ept_order = case_when(
    grepl("Ephemeroptera", TAXON_GROUP_NAME) ~ "Ephemeroptera",
    grepl("Plecoptera",    TAXON_GROUP_NAME) ~ "Plecoptera",
    grepl("Trichoptera",   TAXON_GROUP_NAME) ~ "Trichoptera"
  ))

ept_keys <- ept_taxa$TAXON_LIST_ITEM_KEY

result <- result %>%
  mutate(is_ept = TAXON_LIST_ITEM_KEY %in% ept_keys)
#Check how many taxon records are EPT
table(result$is_ept)

# Map each taxon key to its EPT order
order_lookup <- ept_taxa %>% distinct(TAXON_LIST_ITEM_KEY, ept_order)

# Build per-sample presence/absence for each order
model_data <- result %>%
  left_join(order_lookup, by = "TAXON_LIST_ITEM_KEY") %>%
  group_by(SAMPLE_ID) %>%
  summarise(
    SITE_ID      = first(SITE_ID),
    SAMPLE_DATE  = first(SAMPLE_DATE),
    temp_on_day  = first(temp_on_day),
    ALKALINITY   = first(ALKALINITY),
    DISCHARGE    = first(DISCHARGE),
    Ephemeroptera = any(ept_order == "Ephemeroptera", na.rm = TRUE),
    Plecoptera    = any(ept_order == "Plecoptera",    na.rm = TRUE),
    Trichoptera   = any(ept_order == "Trichoptera",   na.rm = TRUE),
    .groups = "drop"
  )
model_data <- model_data %>%
  filter(!is.na(DISCHARGE)) %>%
  mutate(
   # DISCHARGE = as.factor(DISCHARGE),
    temp_z    = as.numeric(scale(temp_on_day)),
    alk_z     = as.numeric(scale(ALKALINITY))
  )
#Check these are mayfly/stonefly/caddisfly families
result %>% filter(is_ept) %>% distinct(TAXON_NAME) %>% head(20)

#Extra: check that non-EPT looks right too
result %>% filter(!is_ept) %>% distinct(TAXON_NAME) %>% head(20)

#for table summaries
#keep the columns that do NOT contain "_seg"
keep_cols <- grep("_seg", names(aspt_model_data), invert = TRUE, value = TRUE)

#New data.table with only those columns
aspt_model_data_sum <- aspt_model_data[, ..keep_cols]
aspt_model_sum<- as.data.table(summary(aspt_model_data_sum))

#add in the year
aspt_model_data <- aspt_model_data %>%
  mutate(SAMPLE_YEAR = as.integer(format(SAMPLE_DATE, "%Y")))

model_data <- model_data %>%
  mutate(SAMPLE_YEAR = as.integer(format(SAMPLE_DATE, "%Y")))
model_data$Ephemeroptera<- as.numeric(model_data$Ephemeroptera)
model_data$Plecoptera<- as.numeric(model_data$Plecoptera)
model_data$Trichoptera<- as.numeric(model_data$Trichoptera)
##Finally, create the models

#First model: WHPT score response - average score per taxon as the variable
scoremodel<- lmer(WHPT_ASPT ~ temp_z + alk_z +(1|SITE_ID)+ (1|SAMPLE_YEAR) + SAMPLE_YEAR +
               DISCHARGE + temp_z*alk_z, data=aspt_model_data)
#Model of Ephemeroptera presence (mayflies)
ephmodel<- glmer(Ephemeroptera ~ temp_z +alk_z +(1|SITE_ID)+DISCHARGE + (1|SAMPLE_YEAR)  +
                   temp_z*alk_z, data = model_data, family=binomial)
#Model for Trichoptera presence (caddisflies)
trimodel<- glmer(Trichoptera ~ temp_z +alk_z +(1|SITE_ID)+DISCHARGE + (1|SAMPLE_YEAR)  +
                   temp_z*alk_z, data = model_data, family=binomial)
#Model for Plecoptera presence (stoneflies)
plemodel<- glmer(Plecoptera ~ temp_z +alk_z +(1|SITE_ID)+DISCHARGE + (1|SAMPLE_YEAR) +
                   temp_z*alk_z, data = model_data, family=binomial)

#Print the model outputs
summary(scoremodel)
summary(ephmodel)
summary(trimodel)
summary(plemodel)

#Residual checks
#Residuals vs fitted (checks for patterns / non-linearity)
plot(fitted(scoremodel), resid(scoremodel),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, lty = 2, col = "red")

# Q-Q plot (checks normality of residuals)
qqnorm(resid(scoremodel)); qqline(resid(scoremodel))

##Residuals for the binomial models
#Ephemeroptera
simE <- simulateResiduals(fittedModel = ephmodel)

plot(simE)              # two diagnostic plots (QQ + residual vs predicted)
testDispersion(simE)    # formal dispersion test with p-value
testZeroInflation(simE)      # presence/absence data often has zero-inflation

#Trichoptera
simT <- simulateResiduals(fittedModel = trimodel)

plot(simT)              # two diagnostic plots (QQ + residual vs predicted)
testDispersion(simT)    # formal dispersion test with p-value
testZeroInflation(simT)      # presence/absence data often has zero-inflation

#Plecoptera
simP <- simulateResiduals(fittedModel = plemodel)

plot(simP)              # two diagnostic plots (QQ + residual vs predicted)
testDispersion(simP)    # formal dispersion test with p-value
testZeroInflation(simP)      # presence/absence data often has zero-inflation

#####Plots
##Score model plots
#temp
tempmod <- ggpredict(scoremodel, terms = "temp_z")

tempmod$x_real <- tempmod$x * sd(aspt_model_data$temp_on_day, na.rm = TRUE) +
  mean(aspt_model_data$temp_on_day, na.rm = TRUE)

tempmodplot<- ggplot(tempmod, aes(x = x_real, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line() +
  labs(x = "Temperature (°C)", y = "Predicted ASPT") +
  theme_minimal()
tempmodplot

#Alk
alkmod <- ggpredict(scoremodel, terms = "alk_z")

alkmod$x_real <- alkmod$x * sd(aspt_model_data$ALKALINITY, na.rm = TRUE) +
  mean(aspt_model_data$ALKALINITY, na.rm = TRUE)

alkmodplot<- ggplot(alkmod, aes(x = x_real, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line() +
  labs(x = "Alkalinity", y = "Predicted ASPT") +
  theme_minimal()
alkmodplot

#discharge
dismod <- ggpredict(scoremodel, terms = "DISCHARGE")

dismodplot <- ggplot(dismod, aes(x = x, y = predicted)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  labs(x = "Discharge", y = "Predicted ASPT") +
  theme_minimal()
dismodplot

##Ephemeroptera model plots
#Temp
Emod<- ggpredict(ephmodel, terms = 'temp_z')
Emod$x_real <- Emod$x * sd(model_data$temp_on_day, na.rm = TRUE) +
  mean(model_data$temp_on_day, na.rm = TRUE)

Emodplot<- ggplot(Emod, aes(x = x_real, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line() +
  labs(x = "Temperature (°C)", y = "Ephemeroptera") +
  theme_minimal()
Emodplot

#Alk
EmodA<- ggpredict(ephmodel, terms = 'alk_z')
EmodA$x_real <- EmodA$x * sd(model_data$ALKALINITY, na.rm = TRUE) +
  mean(model_data$ALKALINITY, na.rm = TRUE)

EmodAplot<- ggplot(EmodA, aes(x = x_real, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line() +
  labs(x = "Alkalinity", y = "Ephemeroptera") +
  theme_minimal()
EmodAplot

#discharge
EmodD <- ggpredict(ephmodel, terms = "DISCHARGE")

EmodDplot<- ggplot(EmodD, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line() +
  labs(x = "Discharge", y = "Ephemeroptera") +
  theme_minimal()
EmodDplot

##Trichoptera model plot
Tmod<- ggpredict(trimodel, terms = 'temp_z')
Tmod$x_real <- Tmod$x * sd(model_data$temp_on_day, na.rm = TRUE) +
  mean(model_data$temp_on_day, na.rm = TRUE)

Tmodplot<- ggplot(Tmod, aes(x = x_real, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line() +
  labs(x = "Temperature (°C)", y = "Trichoptera") +
  theme_minimal()
Tmodplot

#Alkalinity
TmodA<- ggpredict(trimodel, terms = 'alk_z')
TmodA$x_real <- TmodA$x * sd(model_data$ALKALINITY, na.rm = TRUE) +
  mean(model_data$ALKALINITY, na.rm = TRUE)

TmodAplot<- ggplot(TmodA, aes(x = x_real, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line() +
  labs(x = "Alkalinity", y = "Trichoptera") +
  theme_minimal()
TmodAplot

#Discharge
TmodD <- ggpredict(trimodel, terms = "DISCHARGE")

TmodDplot<- ggplot(TmodD, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line() +
  labs(x = "Discharge", y = "Trichoptera") +
  theme_minimal()
TmodDplot
##Plecoptera model plots
#Temp
Pmod<- ggpredict(plemodel, terms = 'temp_z')
Pmod$x_real <- Pmod$x * sd(model_data$temp_on_day, na.rm = TRUE) +
  mean(model_data$temp_on_day, na.rm = TRUE)

Pmodplot<- ggplot(Pmod, aes(x = x_real, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line() +
  labs(x = "Temperature (°C)", y = "Plecoptera") +
  theme_minimal()
Pmodplot

#Alkalinity
PmodA<- ggpredict(plemodel, terms = 'alk_z')
PmodA$x_real <- PmodA$x * sd(model_data$ALKALINITY, na.rm = TRUE) +
  mean(model_data$ALKALINITY, na.rm = TRUE)

PmodAplot<- ggplot(PmodA, aes(x = x_real, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line() +
  labs(x = "Alkalinity", y = "Plecoptera") +
  theme_minimal()
PmodAplot

#Discharge
PmodD <- ggpredict(plemodel, terms = "DISCHARGE")

PmodDplot<- ggplot(PmodD, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line() +
  labs(x = "Discharge", y = "Plecoptera") +
  theme_minimal()
PmodDplot

##Combo plots
(Emodplot + Tmodplot + Pmodplot) + plot_annotation(title = "Predicted probability of taxa presence against Temperature",
                                                   tag_levels = "A")

(EmodAplot + TmodAplot + PmodAplot) + plot_annotation(title = "Predicted probability Taxa presence against Alkalinity",
                                                   tag_levels = "A") 

(EmodDplot + TmodDplot + PmodDplot) + plot_annotation(title = "Predicted probability Taxa presence against Discharge",
                                                      tag_levels = "A") 

(tempmodplot + alkmodplot + dismodplot) + plot_annotation(title = "Predicted ASPT against predictors",
                                           tag_levels = "A") 

yearscore<- ggplot(aspt_model_data, aes(x=SAMPLE_YEAR, y=WHPT_ASPT)) +geom_point()+
  labs(x='Sample Year', y='ASPT')
yearscore

tempyear<- ggplot(aspt_model_data, aes(x=SAMPLE_YEAR, y=temp_on_day)) +geom_point()+
  labs(x='Sample Year', y='Temperature on sample day (°C)')
tempyear

ephtemp<- ggplot(model_data, aes(x=as.factor(Ephemeroptera), y=temp_on_day))+ geom_point()+
  labs(x='Ephemeroptera presence', y='Temperature (°C)')
ephtemp

ephalk<- ggplot(model_data, aes(x=as.factor(Ephemeroptera), y=ALKALINITY))+ geom_violin()+
  labs(x='Ephemeroptera presence', y='Alkalinity')
ephalk

tritemp<- ggplot(model_data, aes(x=as.factor(Trichoptera), y=temp_on_day))+ geom_point()+
  labs(x='Trichopter presence', y='Temperature (°C)')
tritemp

pletemp<- ggplot(model_data, aes(x=as.factor(Ephemeroptera), y=temp_on_day))+ geom_point()+
  labs(x='Ephemeroptera presence', y='Temperature (°C)')
pletemp

tempalk<- ggplot(aspt_model_data, aes(x=alk_z, y=temp_z)) + geom_point()
tempalk

bmwpwhpt<- ggplot(aspt_model_data, aes(x=BMWP_ASPT, y=WHPT_ASPT)) + geom_point() +geom_smooth()+
  labs(x='BMWP ASPT', y='WHPT ASPT')
bmwpwhpt

dischscore<- ggplot(aspt_model_data, aes(x=as.factor(DISCHARGE), y=WHPT_ASPT)) + geom_boxplot()+
  geom_hline(yintercept = mean(aspt_model_data$WHPT_ASPT, na.rm = TRUE),
             colour = "blue", linewidth = 1)+
  labs(x='Discharge', y='ASPT')
dischscore

long_data <- model_data %>%
  pivot_longer(cols = c(Ephemeroptera, Trichoptera, Plecoptera),
               names_to = "taxon",
               values_to = "present") %>%
  filter(present == TRUE)          # keep only where the taxon was recorded

grouptaxtemp<- ggplot(long_data, aes(x = taxon, y = temp_on_day)) +
  geom_boxplot()+
  labs(x='Order',y='Temperature (°C)' )
grouptaxtemp

grouptaxalk<- ggplot(long_data, aes(x=taxon, y=ALKALINITY)) + geom_boxplot()+
  labs(x='Order', y='Alkalinity')
grouptaxalk
# data.table — replace mean_temp_col with your actual column name
site_means <- result[, .(site_mean_temp = mean(mean_temp, na.rm = TRUE)),
                     by = SITE_ID]

sitetemp<- ggplot(site_means, aes(x = SITE_ID, y = site_mean_temp)) +
  geom_point() +
  labs(x = "Site", y = "Mean temperature (°C)")
sitetemp

matchcoordplot <- ggplot() +
  geom_segment(data = plot_sample,
               aes(x = LON, y = LAT,
                   xend = FNODE_LON_seg, yend = FNODE_LAT_seg),
               colour = "grey60", alpha = 0.6) +
  geom_point(data = plot_sample,
             aes(x = LON, y = LAT, colour = "River coordinate"), size = 0.8) +
  geom_point(data = plot_sample,
             aes(x = FNODE_LON_seg, y = FNODE_LAT_seg, colour = "Temperature coordinate"),
             size = 0.8) +
  scale_colour_manual(name = "Location",
                      values = c("River coordinate" = "blue", "Temperature coordinate" = "red")) +
  coord_quickmap() +
  theme_minimal()
matchcoordplot

taxaNaspt<- ggplot(aspt_model_data, aes(x=WHPT_N_TAXA, y=WHPT_ASPT)) +geom_point() 
taxaNaspt
