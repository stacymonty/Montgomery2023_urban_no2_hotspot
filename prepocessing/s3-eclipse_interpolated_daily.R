# Goal: interpolated surface for ERL submission

# clean-up:
#       1. separate tuning, fitting, and evaluating

# to do:
#       X force TZ to local time
#       + use all days at once (instead of doing by-day)?
#       IDW:
#       X use loocv to get best idp or buffer for idw
#       X map: what does it look like with power = 1?
#       Kriging:
#       + revisit data cleaning -> then revisit variogram
#               e.g. what's going on with cottage + 113th?
#       + automap by day on training data to select the best variogram
#       + try normalizing e.g. removing daily averages + just doing resids?
#       + what if you normalize by location?
#       + try wind-direction aware kriging
#       RF
#       X LOOCV to choose parameters
#       X map: what do results look like?
#       X fit and evaluate by day
#       SVM
#       X LOOCV to choose parameters
#       X map: what do results look like?
#       X fit and evaluate by data
#       Future
#       + include all days (with day, wind dir, speed as inputs)
#       + how to use directionality between sites as part of modeling?
#       + try with coords instead of with dists
#       + switch to flexible bayesian implementation like an adult
#       NN
#       + better selection of neurons (?)
#       evaluation
#       + add EPA sites and values to the automap


# 0. set-up --------------------------------------------------------------------

# packages
library(tidyr)
library(dplyr)
library(sf)
library(lubridate)
library(data.table)
library(gstat)
library(ggplot2)
library(spdep)
library(purrr)
library(sp)
library(patchwork) # collect different ggplots
library(furrr)
future::plan(multisession, workers = 4)

# potential interpolation options
library(caret) # general ml package
library(ranger) # random forests
library(neuralnet) # neural nets
library(kernlab) # support vector machine


# 1.3 km grid
grid <- read_sf("lur/data/processed/all_data_on_chicago_clipped.shp") %>%
        st_transform(4326) %>%
        # make maxdist relevant
        st_transform(32616)

# first calculate distances between all pairs of points in the data
grid_sf <- st_transform(st_centroid(grid)[,"u"], 32616)

# check spatial extent
ggplot(grid) + geom_sf() + theme_void()

source("commtgs/scripts/functions.R")

# change to okay for v slow things to run lol
slow <- "ok"

# 1. eclipse: process ----------------------------------------------------------

# load raw eclipse data: august and february
if(!"eclipse_stacy_qaqc.csv" %in% list.files("lur/data/processed")){
        
        # pull in correct coordinates for EPA
        epa <- fread("../ProjectEclipse_Chicago/commtgs/data/processed/aqs_clean.csv") %>%
                group_by(site) %>%
                slice(n = 1) %>%
                mutate(DeviceFriendlyName = toupper(case_when(
                        site == "173157" ~ "EPA Springfield Pump Station",
                        site == "1731219" ~ "EPA Kennedy Expressway",
                        site == "17311" ~ "EPA Village Garage",
                        site == "173176" ~ "EPA COM ED Maintenance Bldg",
                        site == "17316005" ~ "EPA Liberty School",
                        site == "17314002" ~ "EPA Cicero",
                        TRUE ~ "TBA"))) %>%
                ungroup() %>%
                dplyr::select(DeviceFriendlyName2 = DeviceFriendlyName,
                              correct_lon = longitude,
                              correct_lat = latitude) %>%
                filter(DeviceFriendlyName2 != "TBA")
        
        # get official JCDecaux coordinates
        jcd <- fread("../ProjectEclipse_Chicago/sampling_setup/data/raw/busstops.csv") %>%
                mutate(DeviceFriendlyName2 = toupper(paste(public_nam, " (", dir, ")", sep = ""))) %>%
                dplyr::select(correct_lon = point_x, correct_lat = point_y, DeviceFriendlyName2) %>% 
                # add in Nicolas's list
                rbind(
                        fread("../ProjectEclipse_Chicago/commtgs/data/raw/eclipse_jcd_locations_chicago_20210615.csv") %>%
                                mutate(DeviceFriendlyName2 = toupper(paste(DeviceFriendlyName, " (", Direction, ")", sep = ""))) %>%
                                dplyr::select(correct_lon = Longitude, correct_lat = Latitude, DeviceFriendlyName2)    
                ) %>%
                # adjust names to be comparable
                mutate(DeviceFriendlyName2 = gsub("STREET |PLACE |BOULEVARD |[.]|(W LEG) |DRIVE |AVENUE ", "", DeviceFriendlyName2)) %>%
                mutate(DeviceFriendlyName2 = gsub("(WB) (WB)", "(WB)", DeviceFriendlyName2)) %>%
                mutate(DeviceFriendlyName2 = gsub("(EB) (EB)", "(EB)", DeviceFriendlyName2)) %>%
                mutate(DeviceFriendlyName2 = gsub("(SB) (SB)", "(SB)", DeviceFriendlyName2)) %>%
                mutate(DeviceFriendlyName2 = gsub("(NB) (NB)", "(NB)", DeviceFriendlyName2)) %>%
                mutate(DeviceFriendlyName2 = gsub("SWB", "SB", DeviceFriendlyName2)) %>%
                mutate(DeviceFriendlyName2 = case_when(
                        grepl("POLK & L", DeviceFriendlyName2) ~ 'POLK & L "ABOVE" (WB)',
                        grepl("ST JOSEPH HOSPITAL", DeviceFriendlyName2) ~ "LAKE SHORE DR & JOSEPH HOSP (SB)",
                        TRUE ~ DeviceFriendlyName2
                )) %>%
                # if in both data sets, only keep open data bus stop location
                group_by(DeviceFriendlyName2) %>%
                slice(n = 1) %>%
                ungroup() %>% 
                rbind(epa)
        
        # get meta data
        eclipse_list <- fread("../ProjectEclipse_Chicago/commtgs/data/raw/eclipse_metadata.csv") %>%
                # remove missing and incorrect
                filter(!grepl("Removed|Stolen|Ecopod|soon|Destroyed|Plant|pre-deployment|Alex", DeviceFriendlyName, ignore.case = TRUE)) %>%
                filter(DeviceFriendlyName!="") %>%
                # fix specific incorrect names
                mutate(DeviceFriendlyName = case_when(
                        grepl("Polk & L", DeviceFriendlyName) ~ 'Polk & L "Above" (WB)',
                        grepl("39th & Damen", DeviceFriendlyName) ~ "39th & Damen (WB)",
                        grepl("Sheridan & Belle Plaine", DeviceFriendlyName) ~ "Sheridan & Belle Plaine (NB)",
                        grepl("Vincennes 7 75th|Vincennes 75th", DeviceFriendlyName) ~ "Vincennes & 75th (NB)",
                        grepl("Racine & 78th", DeviceFriendlyName) ~ "Racine & 78th (NB)",
                        TRUE ~ DeviceFriendlyName
                )) %>%
                # link with corrected bus shelter coordinates
                mutate(DeviceFriendlyName2 = toupper(DeviceFriendlyName)) %>%
                mutate(DeviceFriendlyName2 = gsub(" ST | 13 E | DRIVE | AVE | AVENUE ", " ", DeviceFriendlyName2)) %>%
                mutate(DeviceFriendlyName2 = gsub("39TH", "PERSHING", DeviceFriendlyName2)) %>%
                mutate(DeviceFriendlyName2 = gsub("SWB", "SB", DeviceFriendlyName2)) %>%
                # drop EPA endings
                mutate(DeviceFriendlyName2 = ifelse(grepl("EPA", DeviceFriendlyName2),
                                                    gsub(" (A|B|C|D|E|F)$", "", DeviceFriendlyName2), 
                                                    DeviceFriendlyName2)) %>%
                # to do: fix GPS coords for EPA stations
                left_join(jcd) %>%
                # fix missing misc annotations
                mutate(MiscAnnotation = case_when(
                        grepl("EPA", DeviceFriendlyName) ~ "EPA Co-Location",
                        # ones where it didn't get filled in are uniformly SRS
                        DeviceFriendlyName %in% c("Kedzie & 53rd (SB)", "47th & Halsted (WB)",
                                                  "79th St & Aberdeen (EB)", "Halsted & 87th (SB)",
                                                  "Racine & 78th (WB)") ~ "Community Selected",
                        DeviceFriendlyName %in% c("Jackson & Wabash (EB)", "State & 83rd (NB)",
                                                  "Vincennes & 75th (NB)", "Stony Island & 89th St (NB)") ~ "Stratified Random Sample",
                        TRUE ~ MiscAnnotation
                )) %>%
                # fix missing obs
                mutate(correct_lon = ifelse(is.na(correct_lon), Longitude, correct_lon),
                       correct_lat = ifelse(is.na(correct_lat), Latitude, correct_lat)) %>%
                dplyr::select(DeviceFriendlyName, MiscAnnotation, correct_lon, correct_lat) %>%
                # only keep one row for each unique location
                group_by(DeviceFriendlyName) %>%
                slice(n = 1) %>%
                ungroup()
        
        # before slicing:
        # 1. check that you've actually got unique misc annotations for every site
        # test <- eclipse_list %>% group_by(DeviceFriendlyName) %>% mutate(n = length(unique(MiscAnnotation))) %>% arrange(n, DeviceFriendlyName)
        # 2. corrected GPS coords for everything except 136th and o (which only shows up once and looks correct)
        # eclipse_list$DeviceFriendlyName2[!eclipse_list$DeviceFriendlyName2 %in% jcd$DeviceFriendlyName2] %>% unique()
        # 3. check that the locations are reasonable on a map
        # map <- sf::st_as_sf(eclipse_list, coords = c("correct_lon", "correct_lat"), crs = 4326)
        # ggplot(map) + geom_sf()
        
        remove(jcd)
        remove(epa)
        
        # load
        eclipse_raw <- fread("../ProjectEclipse_Chicago/commtgs/data/raw/eclipse.csv")
        
        # fix hourly data (note: this is pretty slow)
        eclipse_raw[, datetime_local := ymd_hms(ReadingDateTimeLocal),]
        eclipse_raw[, date_hour := floor_date(datetime_local, unit = "hours")]
        eclipse_raw[, date_day := floor_date(datetime_local, unit = "days")]
        eclipse_raw[, month := month(date_day)]
        
        # only keep necessary columns
        cols <- c("DeviceFriendlyName", "MSRDeviceNbr", "date_hour", "date_day", "month",
                  "CalibratedPM25", "CalibratedO3", "CalibratedNO2", "CalibrationVersion")
        eclipse_raw <- eclipse_raw[, ..cols]
        
        # filter to first year of observation, devicefriendlynames in eclipse_list
        # this gets rid of all the cases where there's no calibration
        eclipse_raw <- eclipse_raw[!grepl("Removed|Stolen|Ecopod|soon|Destroyed|Plant|pre-deployment|Alex", DeviceFriendlyName, ignore.case = TRUE)]
        eclipse_raw <- eclipse_raw[DeviceFriendlyName!=""]
        eclipse_raw <- eclipse_raw[month %in% c(8, 2)]
        
        # N readings, N unique locations, N sensor days
        table(eclipse_raw$month)
        length(unique(eclipse_raw$DeviceFriendlyName[eclipse_raw$month == 8])) # 123
        length(unique(eclipse_raw$DeviceFriendlyName[eclipse_raw$month == 2])) # 106
        nrow(distinct(eclipse_raw[eclipse_raw$month == 8,c("DeviceFriendlyName", "date_hour")])) # 79,390
        nrow(distinct(eclipse_raw[eclipse_raw$month == 8,c("DeviceFriendlyName", "date_day")])) # 3,363
        nrow(distinct(eclipse_raw[eclipse_raw$month == 2,c("DeviceFriendlyName", "date_hour")])) # 63,494
        nrow(distinct(eclipse_raw[eclipse_raw$month == 2,c("DeviceFriendlyName", "date_day")])) # 2,702
        
        # collapse to hourly averages
        eclipse_raw <- eclipse_raw[, .(CalibratedPM25 = mean(CalibratedPM25),
                                       CalibratedO3 = mean(CalibratedO3),
                                       CalibratedNO2 = mean(CalibratedNO2),
                                       n_hour = .N),
                                   by = c("DeviceFriendlyName", "MSRDeviceNbr", "date_day", "date_hour", "month")]
        
        # keep only if hour is 75% complete
        eclipse_raw <- eclipse_raw[n_hour >= 9]
        
        # N readings, N unique locations, N sensor days
        table(eclipse_raw$month) # 77,256 in August and 61,230 in February
        length(unique(eclipse_raw$DeviceFriendlyName[eclipse_raw$month == 8])) # 121
        length(unique(eclipse_raw$DeviceFriendlyName[eclipse_raw$month == 2])) # 106
        nrow(distinct(eclipse_raw[eclipse_raw$month == 8,c("DeviceFriendlyName", "date_day")])) # 3,311
        nrow(distinct(eclipse_raw[eclipse_raw$month == 2,c("DeviceFriendlyName", "date_day")])) # 2,658
        
        # keep only if day is 75% complete
        eclipse_raw[, n_day:= .N, by = c("DeviceFriendlyName", "MSRDeviceNbr", "date_day")]
        eclipse_raw <- eclipse_raw[n_day >= 18]
        
        # N readings, N unique locations, N sensor days
        table(eclipse_raw$month) # 76,302 in August and 60,021 in February
        length(unique(eclipse_raw$DeviceFriendlyName[eclipse_raw$month == 8])) # 121
        length(unique(eclipse_raw$DeviceFriendlyName[eclipse_raw$month == 2])) # 106
        nrow(distinct(eclipse_raw[eclipse_raw$month == 8,c("DeviceFriendlyName", "date_day")])) # 3,311
        nrow(distinct(eclipse_raw[eclipse_raw$month == 2,c("DeviceFriendlyName", "date_day")])) # 2,658
        
        # fwrite(eclipse_raw, "../ProjectEclipse_Chicago/lur/data/processed/eclipse_hrly_stacy_clean.csv")
        # eclipse_raw <- fread("lur/data/processed/eclipse_hrly_stacy_clean.csv")
        
        # clean up DeviceFriendlyName
        eclipse_raw <- eclipse_raw %>%
                # fix specific incorrect names (N = 131 unique locations)
                mutate(DeviceFriendlyName = case_when(
                        grepl("Polk & L", DeviceFriendlyName) ~ 'Polk & L "Above" (WB)',
                        grepl("39th & Damen", DeviceFriendlyName) ~ "39th & Damen (WB)",
                        grepl("Sheridan & Belle Plaine", DeviceFriendlyName) ~ "Sheridan & Belle Plaine (NB)",
                        grepl("Vincennes 7 75th|Vincennes 75th", DeviceFriendlyName) ~ "Vincennes & 75th (NB)",
                        grepl("Racine & 78th", DeviceFriendlyName) ~ "Racine & 78th (NB)",
                        TRUE ~ DeviceFriendlyName
                )) %>%
                # link with correct lat / lon
                left_join(eclipse_list)
        
        # and get nearest neighbors
        sensors <- distinct(eclipse_raw[,c("DeviceFriendlyName", "correct_lat", "correct_lon")]) %>%
                st_as_sf(coords = c("correct_lon", "correct_lat"), crs = 4326) %>%
                st_transform(32616)
        
        # identify 5 (?) nearest neighbors
        n5 <- knearneigh(sensors, k = 5, use_kd_tree = FALSE)
        n5$nn <- cbind(1:nrow(sensors), n5$nn)
        d5 <- dnearneigh(sensors, d1 = 0, d2 = 5000, use_kd_tree = FALSE)
        
        # for each device and its five nearest neighbors, calculate correlation
        index <- 1:nrow(sensors)
        cors <- map_df(index, check_nbrs_cor)
        
        # what should the cutoff be?
        cor_check <- eclipse_raw %>% 
                group_by(DeviceFriendlyName, month) %>% 
                summarise(n_days = length(unique(date_day))) %>% 
                ungroup() %>% spread(key = month, value = n_days) %>% 
                left_join(cors)
        
        # 0.25 seems like the most defensible based on the distribution
        ggplot(cor_check) + geom_histogram(aes(x = cor_aug)) + geom_histogram(aes(x = cor_feb), fill = "blue", alpha = 0.5) + theme_minimal()
        
        eclipse_raw <- eclipse_raw %>%
                group_by(DeviceFriendlyName, month) %>%
                mutate(n_month = length(unique(date_day))) %>%
                ungroup() %>%
                filter(n_month >= 5)
        
        table(eclipse_raw$month) # 59,883 in August and 75,908 in February
        length(unique(eclipse_raw$DeviceFriendlyName[eclipse_raw$month == 8])) # 104
        length(unique(eclipse_raw$DeviceFriendlyName[eclipse_raw$month == 2])) # 103
        nrow(distinct(eclipse_raw[eclipse_raw$month == 8,c("DeviceFriendlyName", "date_day")])) # 3,205
        nrow(distinct(eclipse_raw[eclipse_raw$month == 2,c("DeviceFriendlyName", "date_day")])) # 2,5355
        
        # process cors
        cors <- cors %>% 
                dplyr::select(-n) %>% 
                gather(key = month, value = cor, cor_aug, cor_feb) %>% 
                mutate(month = as.numeric(ifelse(month == "cor_feb", 2, 8)))
        
        # add in cors
        eclipse_raw <- eclipse_raw %>%
                left_join(cors) %>%
                filter((cor > 0.5 & month == 2) | (cor > 0.5 & month == 8))
        
        table(eclipse_raw$month) # 53,972 in August and 75,147 in February
        length(unique(eclipse_raw$DeviceFriendlyName[eclipse_raw$month == 8])) # 103
        length(unique(eclipse_raw$DeviceFriendlyName[eclipse_raw$month == 2])) # 93
        nrow(distinct(eclipse_raw[eclipse_raw$month == 8,c("DeviceFriendlyName", "date_day")])) # 3,173
        nrow(distinct(eclipse_raw[eclipse_raw$month == 2,c("DeviceFriendlyName", "date_day")])) # 2,286
        
        # export final, clean, hourly data
        fwrite(eclipse_raw, "lur/data/processed/eclipse_stacy_qaqc.csv")
        rm(eclipse_raw)
        eclipse_clean <- fread("lur/data/processed/eclipse_stacy_qaqc.csv") %>%
                mutate(date_hour = force_tz(date_hour, tzone = "America/Chicago"),
                       date_day = force_tz(date_day, tzone = "America/Chicago"))
}else{
        
        eclipse_clean <- fread("lur/data/processed/eclipse_stacy_qaqc.csv") %>%
                mutate(date_hour = force_tz(date_hour, tzone = "America/Chicago"),
                       date_day = force_tz(date_day, tzone = "America/Chicago"))
        
        # 2/6/2023: double check Ns
        # august: 103 locations and 3173 device days
        length(unique(eclipse_clean$DeviceFriendlyName[eclipse_clean$month==8]))
        nrow(distinct(eclipse_clean[eclipse_clean$month==8, c("DeviceFriendlyName", "date_day")]))
        # February: 93 locations and 2,286 device-days
        length(unique(eclipse_clean$DeviceFriendlyName[eclipse_clean$month==2]))
        nrow(distinct(eclipse_clean[eclipse_clean$month==2, c("DeviceFriendlyName", "date_day")]))
}

# load EPA data
aqs_clean <- fread("commtgs/data/processed/aqs_clean.csv") %>%
        mutate(date_hour = force_tz(date_hour, tzone = "America/Chicago"),
               date_day = force_tz(date_day, tzone = "America/Chicago")) %>%
        mutate(month = month(date_day)) %>%
        filter(month %in% c(2, 8) & parameter == "NO2") %>%
        st_as_sf(coords = c("longitude", "latitude"), 
                 crs = 4326, remove = FALSE)


# 2. eclipse: set up validation ------------------------------------------------

# set up validation: EPA AQS sites
# com ed was not operating through august, apparently...
aqs_daily <- aqs_clean %>%
        group_by(site_name, site, latitude, longitude,
                 date_day) %>%
        summarise(epa_no2 = mean(sample_measurement)) %>%
        st_as_sf(coords = c("longitude", "latitude"),
                 crs = 4326) %>%
        st_transform(32616) %>%
        st_intersection(grid) %>%
        dplyr::select(site_name:u)
        
# set up test set: eclipse daily data
eclipse_daily <- eclipse_clean %>%
        # exclude EPA co-locs
        filter(!grepl("EPA", DeviceFriendlyName)) %>%
        group_by(DeviceFriendlyName, correct_lon, correct_lat, MiscAnnotation, 
                 month, date_day) %>%
        summarise(CalibratedNO2 = mean(CalibratedNO2)) %>%
        ungroup() %>%
        st_as_sf(coords = c("correct_lon", "correct_lat"),
                 crs = 4326) %>%
        st_transform(32616) %>%
        # this is super slow for some reason? and not necessary.
        st_intersection(grid) %>%
        dplyr::select(DeviceFriendlyName:u)

rm(eclipse_clean); rm(aqs_clean)

#### delete? -------------------------------------------------------------------
# check for stacy (2/6/2023) - what sites are included?
tot_days <- eclipse_daily %>% 
        group_by(DeviceFriendlyName, month) %>%
        summarise(n_days = length(unique(date_day)),
                  CalibratedNO2 = mean(CalibratedNO2)) %>%
        ungroup() %>%
        mutate(per_days = n_days / max(n_days)) %>%
        st_transform(crs = 4326)

# EPA
tot_days_epa <- aqs_clean %>%
        group_by(site_name, month, longitude, latitude) %>%
        summarise(n_days = length(unique(date_day)),
                  CalibratedNO2 = mean(sample_measurement)) %>%
        ungroup() %>%
        mutate(per_days = n_days / max(n_days)) 

# grid
grid_rf_20230207 <- st_read("lur/data/processed/all_data_grid_eclipse_rf.geojson") %>%
        st_transform(crs = 4326)

tot_days$longitude <- st_coordinates(tot_days)[,1]
tot_days$latitude <- st_coordinates(tot_days)[,2]

grid_rf_cleaner <- grid_rf_20230207 %>%
        dplyr::select(u, CalibratedNO2 = rf_month_august) %>%
        mutate(month = 8) %>%
        rbind(grid_rf_20230207 %>%
                      dplyr::select(u, CalibratedNO2 = rf_month_february) %>%
                      mutate(month = 2))

ggplot() +
        geom_sf(data = grid_rf_cleaner) +
        geom_sf(data = tot_days,
                aes(fill = CalibratedNO2, alpha = per_days), pch = 24) + 
        geom_sf(data = tot_days_epa,
                aes(fill = CalibratedNO2, alpha = per_days), pch = 21) +
        scale_fill_gradient2("",low = 'blue', mid = 'white', high = 'red', midpoint = 9) +
        facet_wrap(~month)

library(gridExtra)

g1 <- ggplot() +
        geom_sf(data = grid_rf_cleaner[grid_rf_cleaner$month==8,],
                aes(fill = CalibratedNO2)) +
        geom_sf(data = tot_days[tot_days$month==8,],
                aes(fill = CalibratedNO2, alpha = per_days), pch = 24) + 
        geom_sf(data = tot_days_epa[tot_days_epa$month==8,],
                aes(fill = CalibratedNO2, alpha = per_days), pch = 21) +
        scale_fill_gradient2("",low = 'blue', mid = 'white', high = 'red', midpoint = 9) +
        facet_wrap(~month) + xlab("August")

g2 <-  ggplot() +
        geom_sf(data = grid_rf_cleaner[grid_rf_cleaner$month==2,],
                aes(fill = CalibratedNO2)) +
        geom_sf(data = tot_days[tot_days$month==2,],
                aes(fill = CalibratedNO2, alpha = per_days), pch = 24) + 
        geom_sf(data = tot_days_epa[tot_days_epa$month==2,],
                aes(fill = CalibratedNO2, alpha = per_days), pch = 21) +
        scale_fill_gradient2("",low = 'blue', mid = 'white', high = 'red', midpoint = 14) +
        facet_wrap(~month) + xlab("February")

grid.arrange(g2, g1, nrow = 1)


# interactive version
library(leaflet)
leaflet() %>% 
        # AQS sites
        addProviderTiles(providers$CartoDB.Positron) %>% 
        addCircles(data = tot_days[tot_days$month==2,],
                   ~longitude, ~latitude, 
                   weight = 8,
                   label = ~per_days, 
                   color = "green") %>%
        addCircles(data = tot_days_epa[tot_days_epa$month==2,],
                   ~longitude, ~latitude, 
                   weight = 8,
                   label = ~per_days, 
                   color = "blue")

leaflet() %>% 
        # AQS sites
        addProviderTiles(providers$CartoDB.Positron) %>% 
        addCircles(data = tot_days[tot_days$month==2,],
                   ~longitude, ~latitude, 
                   weight = 8,
                   label = ~per_days, 
                   color = "green")

#### ---------------------------------------------------------------------------
        


eclipse_daily_aug <- eclipse_daily %>% filter(month == 8 & year(date_day) == 2021)
eclipse_daily_feb <- eclipse_daily %>% filter(month == 2)

# idea: split eclipse by location (80/20)
# better idea would probably be LOOCV...
set.seed(3620)
folds_aug <- groupKFold(eclipse_daily_aug$DeviceFriendlyName, k = 5)
set.seed(3620)
folds_feb <- groupKFold(eclipse_daily_feb$DeviceFriendlyName, k = 5)
#lapply(folds_aug, function(x, y) table(y[x]), y = eclipse_daily_aug$DeviceFriendlyName)


# try fixing feb
#low <- c("79th & Kostner (EB)", "Madison & Loomis (WB)")#, "Ave O & 114th (NB)")

# and create resulting data sets
train <- rbind(eclipse_daily_aug[folds_aug$Fold1,], 
               eclipse_daily_feb[folds_feb$Fold5,] #%>% filter(!(month == 2 & DeviceFriendlyName %in% low))
               )
test <- rbind(eclipse_daily_aug[-folds_aug$Fold1,], 
              eclipse_daily_feb[-folds_feb$Fold5,] #%>% filter(!(month == 2 & DeviceFriendlyName %in% low))
              )

# double check that theres no overlap
testday <- unique(train$date_day)[1]
testdays <- c(unique(train$date_day)[1], unique(train$date_day)[40])
ggplot() + 
        geom_sf(data = grid) +
        geom_sf(data = train[train$date_day %in% testdays,], 
                pch = 21, size = 3, aes(fill = CalibratedNO2)) +
        geom_sf(data = test[test$date_day %in% testdays,], 
                pch = 22, size = 3, aes(fill = CalibratedNO2)) +
        scale_fill_distiller(palette = "Spectral") + theme_void() + 
        facet_wrap(~month)

# 3. eclipse: compare interpolations -------------------------------------------

# 3.0 daily mean imputation
mean_daily <- train %>%
        group_by(date_day) %>%
        summarise(citywide_mean = mean(CalibratedNO2))

# 3.1 IDW
# this is a lot, so comment out for now
# first, use LOOCV to find "best" IDP (lowest RMSE)
# in august, 1 and 1.5 often give basically the same performance
# in february, 1 is much better...?
# 2 or more is almost always worse across months
idp <- 1
if(is.na(idp)){
        eval_idw_loocv <- function(date, params = c(0.5, 0.75, 1, 1.25, 1.5, 2)){
                day_idw <- data.frame(date = date, 
                                      n = 1:nrow(train[train$date_day==date,]), 
                                      power05 = NA, power75 = NA, power1 = NA, power125 = NA, power15 = NA, power2 = NA)
                
                
                for(i in 1:length(params)){
                        day_idw[,i+2] <- map_df(1:nrow(train[train$date_day==date,]),
                                                fit_idw_loocv, power = params[i], day = date)
                }
                
                return(day_idw)
                
        }
        
        # there has to be a faster way to do this...
        tictoc::tic()
        power_idw <- map_df(unique(train$date_day), eval_idw_loocv)
        tictoc::toc()
        
        # check results
        rmse <- function(x){sqrt(mean(x^2))}
        power_tune <- power_idw %>% 
                group_by(date) %>% 
                summarise(power05 = rmse(power05),
                          power75 = rmse(power75),
                          power1 = rmse(power1),
                          power125 = rmse(power125), 
                          power15 = rmse(power15),
                          power2 = rmse(power2)) %>%
                mutate(month = month(date))
        
        # consistently lowest when IDP = 1
        ggplot(power_tune) + 
                geom_line(aes(x = date, y = power05)) +
                geom_line(aes(x = date, y = power75), col = "purple") +
                geom_line(aes(x = date, y = power1), col = "blue") + 
                geom_line(aes(x = date, y = power125), col = "red") + 
                geom_line(aes(x = date, y = power15), col = "green") + 
                facet_wrap(~month, scales = "free")
}

# fit idws by day (slow!)
# using future_map_dfr brings it down to one minute :)
if(slow == "ok"){
        tictoc::tic()
        days <- unique(train$date_day)
        idw_daily <- future_map_dfr(days, 
                                    fit_idw_cv,
                                    variable = "CalibratedNO2",
                                    power = 1,
                                    .options = furrr_options(seed = TRUE),
                                    .progress = TRUE) 
        tictoc::toc()
        # takes five minutes
        # consider exporting this so you don't always have a long task
}

#### 3.2 ordinary kriging ------------------------------------------------------

# not worth the brain power right now
# We also investigated ordinary kriging, a geostatistical method commonly used for the interpolation of pollution data (cite).
# However, the approach performed particularly poorly, failing to outperform the citywide mean, and thus we do not report results here.

# to do: just use automap

#### 3.3 Random Forests --------------------------------------------------------

fit_rf_loocv <- function(day = testday, variable = "CalibratedNO2"){
        
        # just pull data for the relevant day
        dat <- train[train$date_day==day,] %>%
                dplyr::select(DeviceFriendlyName, date_day, CalibratedNO2)
        
        # write a fxn to return clean input for ml fxns
        ai_input <- get_input(dat)
        
        # set up tuning grid
        trctrl <- trainControl(method = "LOOCV")
        # parameters are mtry, splitrule, and min.node.size
        params_rf <- expand.grid(#num.trees = 1000?
                                        mtry = c(1, 2, 3, 5, 10),
                                        min.node.size = c(5, 10, 25, 50),
                                        splitrule = "variance"
                                        )
        set.seed(3620)
        tictoc::tic()
        pred_rf <- train(
                z ~ .,
                data = data.frame(ai_input[1]),
                method = "ranger",
                tuneGrid = params_rf,
                trControl = trctrl
                
        )
        tictoc::toc()
        
        return(data.frame(pred_rf$results, date = day))
        
}

# loop through each day to find the best params (and avg rmse)
# estimated to take about an hour
# using furrr instead of purrr to speed things up (massively helps!)
params_rf <- data.frame(mtry = 3, min.node.size = 10)
if(is.na(params_rf[1])){
        tictoc::tic()
        tune_rf <- future_map_dfr(unique(train$date_day), fit_rf_loocv,
                                  .options = furrr_options(seed = TRUE),
                                  .progress = TRUE)
        tictoc::toc()
        
        # what parameters lead to the best performance
        ggplot(tune_rf) + 
                geom_line(aes(x = mtry, y = RMSE, col = factor(date))) + geom_point(aes(x = mtry, y = RMSE, col = factor(date))) +
                facet_wrap(~min.node.size, nrow = 1) +
                theme(legend.position = "none")
        
        # take best result
        tune_rf %>% 
                group_by(min.node.size, mtry) %>% 
                summarise(param = mean(RMSE)) %>% 
                ungroup() %>% 
                # searching within 1SD so we don't end up with something bonkers
                filter(param < min(param) + sd(param))
}


# write fxn to fit rf by day
fit_rf <- function(day = testday, params = parms_rf){
        
        # just pull data for the relevant day
        dat <- train[train$date_day==day,] %>%
                dplyr::select(DeviceFriendlyName, date_day, CalibratedNO2)
        
        # write a fxn to return clean input for ml fxns
        ai_input <- get_input(dat)
        
        # fit the model
        mod_rf <- ranger( z ~ .,
                         data = data.frame(ai_input[1]), 
                         num.trees = 500, 
                         mtry = params_rf$mtry,
                         min.node.size = params_rf$min.node.size)
        
        # predict on other grid cells
        pred_rf <- predict(mod_rf, data.frame(ai_input[2])) %>%
                as_tibble()
        
        # rescale
        pred_rf$z <- denormalize(
                pred_rf$prediction, 
                # this is by daily min and max
                # probably better to use overall values?
                top = max(dat$CalibratedNO2), 
                bottom = min(dat$CalibratedNO2)
        )
        
        # set up for export
        pred_rf <- pred_rf %>%
                mutate(model = "Random Forest") %>%
                dplyr::select(z, model)
        
        pred_rf$u <- grid$u
        
        pred_rf$date <- day
        
        return(pred_rf)
        
}

# test + map one or two 
rf_check <- fit_rf(unique(train$date_day)[1])

grid %>%
        left_join(rf_check) %>%
        ggplot() + geom_sf(aes(fill = z)) +
        scale_fill_distiller(palette = "Spectral")

# fit by day
# lol this takes 7 seconds
# update: restarted RStudio and all the parallelization got bad :(
if(slow == "ok"){
        tictoc::tic()
        rf_daily <- future_map_dfr(unique(train$date_day), 
                            fit_rf,
                            .options = furrr_options(seed = TRUE),
                            .progress = TRUE) 
        tictoc::toc()
}


#### 3.4 SVMs ------------------------------------------------------------------

# main hyperparameter is the kernel
# idea is to pick a transformation such that the observations become more easily (linearly) separable after the 
# 1. linear kernel; has a cost parameter C that has to be tuned
# 2. polynomial kernel: r 0, gamma fixed (e.g. 1/n); tune cost and d - between values 1 and 10
# 3. radial kernel: most used and successful
# cost parameter C (range 2^-5 to 2^15) and gamma in range 2^-15 to 2^3
# mlrHyperoprt R package should help tune

# TO DO
# X use caret to do LOOCV correctly?
# need to choose a kernel, cost fxn, and sigma parameter
# https://www.csie.ntu.edu.tw/~cjlin/papers/guide/guide.pdf
# linear, polynomial, radial basis function, and sigmoid kernels are common
# 1. transform data; scaling; consider RBF; use CV to find best C and sigma; use best to train whole set; test
# need to scale to avoid attributes in greater numeric ranges dominating small ranges - linearly scale to -1, 1
# make sure you scale traing and testing data the same way
# radial basis function kernel can handle nonlinear relationships between class labels and attributes
# (linear kernel is a special case of RBF, and sigmoid kernel behaves similarly)
# polynomial kernel has more hyperparameters, which means more complex model selection

# grid-search: C and y tried in pairs; use the one with the best CV accuracy
# start with exponentially growing sequences of C and sigma e.g. C = 2-5, 2-3, ...2^15; and same for y
# you can do naive grid search
# since there are only two parameters, the time isn't that much
# also should be easy to parallelize

# start with a coarse grid, then hone in on the finer grid
# if you have a really big data set, you can choose a subset / do grid-search / then try with all data

# available models and parameters: http://topepo.github.io/caret/train-models-by-tag.html#support-vector-machines


fit_svm_loocv <- function(day = testday, variable = "CalibratedNO2"){
        
        # just pull data for the relevant day
        dat <- train[train$date_day==day,] %>%
                dplyr::select(DeviceFriendlyName, date_day, CalibratedNO2)
        
        # write a fxn to return clean input for ml fxns
        ai_input <- get_input(dat)
        
        # try with a radial basis function kernel
        trctrl <- trainControl(method = "LOOCV")
        params_rbf <- expand.grid(C = c(2^(-4), 2^(-1), 2, 2^2, 2^4),
                              sigma = c(2^(-5), 2^(-1), 2, 2^2, 2^5))
        set.seed(3620)
        tictoc::tic()
        svm_rbf <- train(
                z ~ .,
                data = data.frame(ai_input[1]),
                method = "svmRadial",
                tuneGrid = params_rbf,
                trControl = trctrl
                
        )
        tictoc::toc()
        
        #data.frame(svm_rbf$results)
        
        tictoc::tic()
        params_poly <- expand.grid(C = c(2^(-5), 2^(-1), 2, 2^5),
                                 degree = c(1, 2, 3),
                                 scale = c(0.001, 0.1, 1, 10))
        svm_poly <- train(
                z ~ .,
                data = data.frame(ai_input[1]),
                method = "svmPoly",
                tuneGrid = params_poly,
                trControl = trctrl
                
        )
        tictoc::toc()
        
        result <- data.frame(svm_poly$results, model = "polynomial svm", date = day, sigma = NA) %>%
                dplyr::select(sigma, C, degree, scale,
                              RMSE:date) %>%
                rbind(data.frame(svm_rbf$results, model = "rbf svm", date = day, degree = NA, scale = NA) %>%
                              dplyr::select(sigma, C, degree, scale,
                                            RMSE:date))
        
        return(result)
}

# for polynomia, maybe C = 0.5, degree = 1, scale = 0.1 (but what does degree 1 even mean)
# very clearly worsens with increasing degree
params_svm <- data.frame(C = 0.5, sigma = 0.0312, model = "rbfdot")
if(is.na(params_svm[1])){
        tictoc::tic()
        tune_svm <- future_map_dfr(unique(train$date_day), fit_svm_loocv,
                                  .options = furrr_options(seed = TRUE),
                                  .progress = TRUE)
        tictoc::toc()
        
        # what parameters lead to the best performance
        ggplot(tune_svm[tune_svm$model=="rbf svm",],
               aes(x = sigma, y = RMSE, group = factor(date))) + 
                geom_point() + geom_line() +
                facet_wrap(~C, nrow = 1) +
                theme(legend.position = "none")
        
        # take best result
        # rbf svm: C = 0.5; sigma = 0.0312
        tune_svm %>% 
                filter(model == "rbf svm") %>%
                group_by(C, sigma) %>% 
                summarise(param = mean(RMSE)) %>% 
                ungroup() %>% 
                # searching within 1SD so we don't end up with something bonkers
                filter(param < min(param) + sd(param))
}

# write fxn to fit svm by day
fit_svm <- function(day = testday, params = params_svm){
        
        # just pull data for the relevant day
        dat <- train[train$date_day==day,] %>%
                dplyr::select(DeviceFriendlyName, date_day, CalibratedNO2)
        
        # write a fxn to return clean input for ml fxns
        ai_input <- get_input(dat)
        
        mod_svm <- ksvm(
                z ~ .,
                data = data.frame(ai_input[1]),
                kernel = params_svm$model, # this is the decision point; i've been using polynomial
                C = params_svm$C, # A parameter to penalize overfitting
                sigma = params_svm$sigma
                #kernel = "polydot",
                #C = 0.3125,
                #scale = 0.1,
                #degree = 2
        )
        
        # predict across the grid
        pred_svm <- predict(mod_svm, data.frame(ai_input[2])) %>%
                as_tibble()
        
        # convert back to NO2 values
        pred_svm$z <- denormalize(
                pred_svm$V1, 
                top = max(dat$CalibratedNO2), 
                bottom = min(dat$CalibratedNO2)
        )
        
        # save grid values
        pred_svm <- pred_svm %>%
                mutate(model = "Support Vector Machine") %>%
                dplyr::select(z, model)
        pred_svm$u <- grid$u
        pred_svm$date <- day
        
        # might consider doing nn and rf in this same fxn to save time
        return(pred_svm)
        
}

# test + map one or two 
# definitely looking wonky
svm_check <- fit_svm(unique(train$date_day)[1])
grid %>%
        left_join(svm_check) %>%
        ggplot() + geom_sf(aes(fill = z)) +
        scale_fill_distiller(palette = "Spectral")

# fit by day
# lol this takes 6 seconds
# update: something broke and now it takes 40 seconds :(
if(slow == "ok"){
        tictoc::tic()
        svm_daily <- future_map_dfr(unique(train$date_day), 
                                   fit_svm,
                                   .options = furrr_options(seed = TRUE),
                                   .progress = TRUE) 
        tictoc::toc()
}

#### 3.5 Neural Nets -----------------------------------------------------------

#fit_nn_loocv <- function(day = testday)
#
#mod_nn <- neuralnet(
#        z ~., # col z is dependent var; all others are ind vars
#        data = nn_input,
#        hidden = c(
#                200, 100, 50, 150
#        ) # four layers of several arbitrarily chosen neurons
#        # what function is getting used in the final transformations
#)
## predict for each grid point
#pred_nn <- predict(mod_nn, grid_dnorm) %>%
#        as_tibble()
## clean tibble of results
#pred_nn$z <- denormalize(
#        pred_nn$V1, 
#        top = max(train_coords$CalibratedNO2), 
#        bottom = min(train_coords$CalibratedNO2)
#)
#pred_nn <- pred_nn %>%
#        mutate(model = "Neural Network") %>%
#        dplyr::select(z, model)
#pred_nn$u <- grid$u

# 4. evaluate results ----------------------------------------------------------

# i think i need this for the ML evals
test_sf <- st_transform(test[test$date_day == testday,],
                        32616)
aqs_sf <- st_transform(aqs_daily[aqs_daily$date_day == testday,],
                       32616)

# set up data: daily estimates by location
test_eval <- test %>%
        # link with citywide averages (by day)
        left_join(st_drop_geometry(mean_daily)) %>%
        # link with idw pred (by u and day)
        left_join(st_drop_geometry(idw_daily) %>%
                          dplyr::select(date_day, u, idw_no2 = var1.pred)) %>%
        # to do: link with kriging
        # link with svm pred (by u and day)
        left_join(svm_daily %>%
                          dplyr::select(date_day = date, u, svm_no2 = z)) %>%
        # link with rf pred (by u and day)
        left_join(rf_daily %>%
                          dplyr::select(date_day = date, u, rf_no2 = z)) 
        # to do: link with neural net output

# aqs: daily estimates by location
aqs_eval_day <- aqs_daily %>%
        # link with citywide averages (by day)
        left_join(st_drop_geometry(mean_daily)) %>%
        # link with idw pred (by u and day)
        left_join(st_drop_geometry(idw_daily) %>%
                          dplyr::select(date_day, u, idw_no2 = var1.pred)) %>%
        # to do: link with kriging
        # link with svm pred (by u and day)
        left_join(svm_daily %>%
                          dplyr::select(date_day = date, u, svm_no2 = z)) %>%
        # link with rf pred (by u and day)
        left_join(rf_daily %>%
                          dplyr::select(date_day = date, u, rf_no2 = z)) 
        # to do: link with neural net output

# aqs overall (month) evaluation
aqs_eval_month <- aqs_eval_day %>%
        # group by site, month
        mutate(month = month(date_day)) %>%
        group_by(site_name, site, month) %>%
        # calculate averages for means, idw, ...
        summarise(epa_no2 = mean(epa_no2),
                  citywide_mean = mean(citywide_mean),
                  idw_no2 = mean(idw_no2),
                  svm_no2 = mean(svm_no2),
                  rf_no2 = mean(rf_no2))

# evaluate results (eclipse CV)
# worth stratifying by month - august is fairly solid; feb is just noise
eval_result(test_eval$CalibratedNO2, test_eval$citywide_mean)
eval_result(test_eval$CalibratedNO2, test_eval$idw_no2)
eval_result(test_eval$CalibratedNO2, test_eval$svm_no2)
eval_result(test_eval$CalibratedNO2, test_eval$rf_no2)

# month averages?
test_avg <- test_eval %>%
        group_by(DeviceFriendlyName, month, u) %>%
        summarise(CalibratedNO2 = mean(CalibratedNO2),
                  count = n(),
                  citywide_mean = mean(citywide_mean),
                  idw_no2 = mean(idw_no2),
                  svm_no2 = mean(svm_no2),
                  rf_no2 = mean(rf_no2)) %>%
        group_by(month) %>%
        # all sensors are in "keep" (N = 15)
        mutate(keep = ifelse(count > 0.75*max(count), 1, 0))

# omg, r2 is defensible! and nrmse is really good!
eval_result(test_avg$CalibratedNO2[test_avg$month==8], test_avg$citywide_mean[test_avg$month==8])
eval_result(test_avg$CalibratedNO2[test_avg$month==8], test_avg$idw_no2[test_avg$month==8])
eval_result(test_avg$CalibratedNO2[test_avg$month==8], test_avg$svm_no2[test_avg$month==8])
eval_result(test_avg$CalibratedNO2[test_avg$month==8], test_avg$rf_no2[test_avg$month==8])

eval_result(test_avg$CalibratedNO2[test_avg$month==2], test_avg$citywide_mean[test_avg$month==2])
eval_result(test_avg$CalibratedNO2[test_avg$month==2], test_avg$idw_no2[test_avg$month==2])
eval_result(test_avg$CalibratedNO2[test_avg$month==2], test_avg$svm_no2[test_avg$month==2])
eval_result(test_avg$CalibratedNO2[test_avg$month==2], test_avg$rf_no2[test_avg$month==2])

eval_result(aqs_eval_month$epa_no2[aqs_eval_month$month==8], aqs_eval_month$citywide_mean[aqs_eval_month$month==8])
eval_result(aqs_eval_month$epa_no2[aqs_eval_month$month==8], aqs_eval_month$idw_no2[aqs_eval_month$month==8])
eval_result(aqs_eval_month$epa_no2[aqs_eval_month$month==8], aqs_eval_month$svm_no2[aqs_eval_month$month==8])
eval_result(aqs_eval_month$epa_no2[aqs_eval_month$month==8], aqs_eval_month$rf_no2[aqs_eval_month$month==8])

eval_result(aqs_eval_month$epa_no2[aqs_eval_month$month==2], aqs_eval_month$citywide_mean[aqs_eval_month$month==2])
eval_result(aqs_eval_month$epa_no2[aqs_eval_month$month==2], aqs_eval_month$idw_no2[aqs_eval_month$month==2])
eval_result(aqs_eval_month$epa_no2[aqs_eval_month$month==2], aqs_eval_month$svm_no2[aqs_eval_month$month==2])
eval_result(aqs_eval_month$epa_no2[aqs_eval_month$month==2], aqs_eval_month$rf_no2[aqs_eval_month$month==2])



# evaluate results (aqs by day) <- idw gives us a small improvement in rmse and nrmse, regardless of site
eval_result(aqs_eval_day$epa_no2, aqs_eval_day$citywide_mean)
eval_result(aqs_eval_day$epa_no2, aqs_eval_day$idw_no2)
eval_result(aqs_eval_day$epa_no2, aqs_eval_day$svm_no2)
eval_result(aqs_eval_day$epa_no2, aqs_eval_day$rf_no2)

for(i in unique(aqs_eval_day$site_name)){
        citywide <- eval_result(aqs_eval_day$epa_no2[aqs_eval_day$site_name==i], 
                                aqs_eval_day$citywide_mean[aqs_eval_day$site_name==i])
        idw <- eval_result(aqs_eval_day$epa_no2[aqs_eval_day$site_name==i], 
                    aqs_eval_day$idw_no2[aqs_eval_day$site_name==i])
        svm <- eval_result(aqs_eval_day$epa_no2[aqs_eval_day$site_name==i], 
                           aqs_eval_day$svm_no2[aqs_eval_day$site_name==i])
        rf <- eval_result(aqs_eval_day$epa_no2[aqs_eval_day$site_name==i], 
                           aqs_eval_day$rf_no2[aqs_eval_day$site_name==i])
        
        print(rbind(citywide, idw, svm, rf))
}

# and evaluate (aqs by month) <- citywide does better at Com Ed, idw does better overall + at Kennedy Expressway
eval_result(aqs_eval_month$epa_no2, aqs_eval_month$citywide_mean)
eval_result(aqs_eval_month$epa_no2, aqs_eval_month$idw_no2)
eval_result(aqs_eval_month$epa_no2, aqs_eval_month$svm_no2)
eval_result(aqs_eval_month$epa_no2, aqs_eval_month$rf_no2)
for(i in unique(aqs_eval_day$site_name)){
        citywide <- eval_result(aqs_eval_month$epa_no2[aqs_eval_month$site_name==i], 
                                aqs_eval_month$citywide_mean[aqs_eval_month$site_name==i])
        idw <- eval_result(aqs_eval_month$epa_no2[aqs_eval_month$site_name==i], 
                           aqs_eval_month$idw_no2[aqs_eval_month$site_name==i])
        svm <- eval_result(aqs_eval_month$epa_no2[aqs_eval_month$site_name==i], 
                           aqs_eval_month$svm_no2[aqs_eval_month$site_name==i])
        rf <- eval_result(aqs_eval_month$epa_no2[aqs_eval_month$site_name==i], 
                           aqs_eval_month$rf_no2[aqs_eval_month$site_name==i])
        print(rbind(citywide, idw, svm, rf))
}
        


# visual check
# lol this seems quite driven by day-to-day variation
ggplot(test_eval, aes(x = CalibratedNO2, y = idw_no2, col = factor(date_day))) + 
        geom_point() + geom_smooth(method = "lm", se = FALSE) + theme(legend.position = "none")

# doing well for august and terribly for feb?
ggplot(test_eval, aes(x = CalibratedNO2, y = rf_no2, col = factor(date_day))) + 
        geom_point() + geom_smooth(method = "lm", se = FALSE) + theme(legend.position = "none") + facet_wrap(~month(date_day))

# monthly averages by grid cell
monthly_interps <- st_drop_geometry(idw_daily) %>%
        mutate(model = "inverse distance weighted") %>%
        dplyr::select(CalibratedNO2 = var1.pred, model, u, date = date_day) %>%
        rbind(svm_daily %>%
                      dplyr::select(CalibratedNO2 = z, model, u, date)) %>%
        rbind(rf_daily %>%
                      dplyr::select(CalibratedNO2 = z, model, u, date)) %>%
        mutate(month = month(date)) %>%
        group_by(month, model, u) %>%
        summarise(CalibratedNO2 = mean(CalibratedNO2))

monthly_points <- train %>%
        mutate(type = "Training Set") %>%
        rbind(test %>% 
                      mutate(type = "Test Set")) %>%
        rbind(aqs_daily %>%
                      mutate(MiscAnnotation = NA,
                             month = month(date_day)) %>%
                      dplyr::select(
                              DeviceFriendlyName = site_name,
                              MiscAnnotation,
                              month, date_day, 
                              CalibratedNO2 = epa_no2,
                              u) %>%
                      mutate(type = "EPA Validation")) %>%
        group_by(DeviceFriendlyName, month, u, type) %>%
        summarise(n_obs = n(),
                  CalibratedNO2 = mean(CalibratedNO2)) %>%
        group_by(month) %>%
        mutate(keep = ifelse(n_obs > 0.75*max(n_obs), 1, 0))

all_interps <- grid %>%
        dplyr::select(u) %>%
        distinct() %>%
        left_join(monthly_interps)
        

# august is defensible, beautiful, good statistics
jpeg("lur/output/hotspots_interpolation_figc1_august.jpg", height = 850, width = 1500, res = 150)
ggplot() + 
        geom_sf(data = all_interps %>% filter(month == 8), 
                aes(fill = CalibratedNO2)) +
        geom_sf(data = monthly_points %>% filter(month == 8 & keep == 1),
                aes(fill = CalibratedNO2, pch = type), size = 3) +
        #scale_fill_distiller(palette = "Spectral") +
        scale_fill_gradient2("",low = 'blue', mid = 'white', high = 'red', midpoint = 9) +
        scale_shape_manual("", values = c(24, 22, 21)) +
        facet_wrap(~ model) + theme_void() +
        theme(legend.position = "bottom")
dev.off()

jpeg("lur/output/hotspots_interpolation_figc2_february.jpg", height = 850, width = 1500, res = 150)
ggplot() + 
        geom_sf(data = all_interps %>% filter(month == 2), 
                aes(fill = CalibratedNO2)) +
        geom_sf(data = monthly_points %>% filter(month == 2 & keep == 1),
                aes(fill = CalibratedNO2, pch = type), size = 3) +
        #scale_fill_distiller(palette = "Spectral") +
        scale_fill_gradient2("",low = 'blue', mid = 'white', high = 'red', midpoint = 14) +
        scale_shape_manual("", values = c(24, 22, 21)) +
        facet_wrap(~ model) + theme_void() +
        theme(legend.position = "bottom")
dev.off()

for(i in seq(1, length(unique(idw_daily$date_day)), 5)){
        g <- ggplot() + geom_sf(data = idw_daily[idw_daily$date_day == unique(idw_daily$date_day)[i],], 
                           aes(fill = var1.pred)) + scale_fill_distiller(palette = "Spectral") + 
                geom_sf(data = train[train$date_day== unique(idw_daily$date_day)[i],],
                        aes(fill = CalibratedNO2), pch  = 21) +
                geom_sf(data = test[test$date_day== unique(idw_daily$date_day)[i],],
                        aes(fill = CalibratedNO2), pch  = 22) +
                geom_sf(data = aqs_daily[aqs_daily$date_day== unique(idw_daily$date_day)[i],],
                        aes(fill = epa_no2), pch  = 24) + theme_minimal() + theme(legend.position = "bottom" )  
        print(g)
}

#### 5. export final results ---------------------------------------------------

fit_rf_final <- function(day = testday, params = parms_rf){
        
        # just pull data for the relevant day
        dat <- eclipse_daily[eclipse_daily$date_day==day,] %>%
                dplyr::select(DeviceFriendlyName, date_day, CalibratedNO2)
        
        # write a fxn to return clean input for ml fxns
        ai_input <- get_input(dat)
        
        # fit the model
        mod_rf <- ranger( z ~ .,
                          data = data.frame(ai_input[1]), 
                          num.trees = 500, 
                          mtry = params_rf$mtry,
                          min.node.size = params_rf$min.node.size)
        
        # predict on other grid cells
        pred_rf <- predict(mod_rf, data.frame(ai_input[2])) %>%
                as_tibble()
        
        # rescale
        pred_rf$z <- denormalize(
                pred_rf$prediction, 
                # this is by daily min and max
                # probably better to use overall values?
                top = max(dat$CalibratedNO2), 
                bottom = min(dat$CalibratedNO2)
        )
        
        # set up for export
        pred_rf <- pred_rf %>%
                mutate(model = "Random Forest") %>%
                dplyr::select(z, model)
        
        pred_rf$u <- grid$u
        
        pred_rf$date <- day
        
        return(pred_rf)
        
}

# 1. run daily rf with complete data set
params_rf <- data.frame(mtry = 3, min.node.size = 10)
if(slow == "ok"){
        tictoc::tic()
        rf_daily_final <- future_map_dfr(unique(train$date_day), 
                                   fit_rf_final,
                                   .options = furrr_options(seed = TRUE),
                                   .progress = TRUE) 
        tictoc::toc()
}

# 2. calculate monthly average values
rf_monthly_final <- rf_daily_final %>%
        mutate(month = ifelse(month(date)==8, "rf_month_august",
                              "rf_month_february")) %>%
        group_by(month, u) %>%
        summarise(CalibratedNO2 = mean(z)) %>%
        ungroup() %>%
        dplyr::select(u, month, CalibratedNO2) %>%
        spread(key = month, value = CalibratedNO2)

# 3. export on the grid
old_grid <- read_sf("lur/data/processed/all_data_on_chicago_clipped.shp")
old_grid %>%
        left_join(rf_monthly_final) %>%
        write_sf("lur/data/processed/all_data_grid_eclipse_rf.geojson")

# 4. and export the daily values
rf_daily_final %>%
        dplyr::select(rf_no2 = z,
                      u,
                      date) %>%
        fwrite("lur/data/processed/daily_data_grid_eclipse_rf.csv")


#### check: with IDW = 1, do we still get meaningful getis ord g*? -------------
# https://geo200cn.github.io/spatialautocorrelation.html#Local_spatial_autocorrelation

# idw values by month
grid_rf <- grid[,c("u")] %>%
        left_join(rf_monthly_final)

## create a nearest neighbors object
grid_nb <- poly2nb(grid_rf, queen = T)
grid_nb_self <- include.self(grid_nb)
grid_w_self <- nb2listw(grid_nb_self, style = "W", zero.policy = TRUE)
#
## try calculating the z score
localgstar_aug <- localG(grid_rf$rf_month_august, grid_w_self, zero.policy = TRUE)
localgstar_feb <- localG(grid_rf$rf_month_february, grid_w_self, zero.policy = TRUE)
#
grid_rf <- mutate(grid_rf,
                    rf_gstar_aug = as.numeric(localgstar_aug),
                    rf_gstar_feb = as.numeric(localgstar_feb))

ggplot() + geom_sf(data = grid_rf, fill = "white") +
        geom_sf(data = grid_rf[grid_rf$rf_gstar_aug >= 1.96,], fill = "red") +
        geom_sf(data = grid_rf[grid_rf$rf_gstar_aug < -1.96,], fill = "blue")

ggplot() + geom_sf(data = grid_rf, fill = "white") +
        geom_sf(data = grid_rf[grid_rf$rf_gstar_feb >= 1.96,], fill = "red") +
        geom_sf(data = grid_rf[grid_rf$rf_gstar_feb < -1.96,], fill = "blue")


# export for stacy
old_grid <- read_sf("lur/data/processed/all_data_on_chicago_clipped.shp")
old_grid %>%
        left_join(st_drop_geometry(grid_idw1[,c("u","idw_month_8", "idw_month_2")])) %>%
        write_sf("lur/data/processed/all_data_grid_eclipse_idw_idp1.geojson")



# notes:
# deterministic interpolation: inverse distance weighting, radial basis functions (multiquadratic and thin plate splines), local polynomial interp
# geostatistical interpolation: ordinary kriging, universial kriging, co-kriging, and regression kriging
# metrics: RMSE, MAE, ME, and corr of real and predicted values

# IDW: linear-weighted combination of points within a range
#       nice for rapid estimation of spatial interpolation (Miras-Avalos et al. 2009)

# radial basis function: "minimizes the total curvature of the surface"
# Zhat(S_o) = sum w_i phi(r_i) + m
# interpolation value at s_o is the RBF based on r_i, radial distance from point s_o to s_i; w_i is the weight; m is the bias value

# local polynomial interpolation: many polynomials within specified overlapping neighborhods
# applicable to surfaces that capture short-range variation (Wang et al. 2014)

# ordinary kriging: stochastic method that assesses spatial correlation between measured points using the semivariance
# assumption: constant unkon mean over the sear neighborhood of x_0
# semivariance is a function over distance that can be isotropic (same by direction) or anisotropic (different)
# cressie 1985: semivariance fxn models different values that can't be realistically measured
# "primarily uses a spherical model when a nugget effect is important, and there is an obvious range and sill effect (Burrough et al. 2015; Cressie 2015)"
# "as a simple prediction tool, OK is considered flexible and relatively accurate (Smith et al. 2005)"

# universal kriging:
# "UK is more appropriate than OK when a geographic trend matches the polynomial function that UK generally uses" (?)
# assumes a general polynomial trend model

# co-kriging, residual kriging
# co-kriging: versatile stat method when there are highly correlated interest and auxiliary variables; uses OK on correlation between y and x
# regression kriging: combination of linear regression and kriging - OLS and then OK with resids
# applicable when residuals are spatially correlated

# wikipedia: "kriging / gaussian process regression: method of interpolation based on gaussian process governed by prior covariances"
# "under suitable assumptions of the prior, kriging gives the best linear unbiased prediction at unsampled locations"
# "kriging can also be understood as a form of bayesian optimization. kriging starts with a prior distribution over functions.
# the prior takes the form of a gaussian process: N samples from a function will be normally distributed

#### AI approach ---------------------------------------------------------------
# https://swilke-geoscience.net/post/2020-09-10-kriging_with_r/kriging/
# https://swilke-geoscience.net/post/spatial_ml/
# "one cannot apply machine learning algorithms to spatial problems out of the box" bc coordinates are autocrrelated
# workaround: instead of describing the data as a function of coords x and y, we focus on distances between coords

day <- ymd("2022-08-15")

# first calculate distances between all pairs of points in the data
train_coords <- train[train$date_day == day,]
grid_sf <- st_transform(st_centroid(grid)[,"u"], 32616)
sample_sf <- st_transform(train_coords, 32616)
test_sf <- st_transform(test[test$date_day == day,],
                        32616)
aqs_sf <- st_transform(aqs_daily[aqs_daily$date_day == day,],
                       32616)

sample_dm <- st_distance(sample_sf, sample_sf) %>%
        as_tibble()
sample_dm <- units::drop_units(sample_dm)

# then calculate distances between all points in the data and all centroids of the prediction grid
grid_dm <- st_distance(grid_sf, sample_sf) %>%
        as_tibble() 
grid_dm <- units::drop_units(grid_dm)

# need to normalize the distances for ML algorithms
normalize <- function(x, bottom, top){
        (x - bottom) / (top - bottom)
}
denormalize <- function(x, bottom, top){
        (top - bottom)*x + bottom
}

# "rather than just normalizing each vector individually"
# "we make an informed decision about the distances that can possibly occur"
# so distances are weighted the same everywhere
bottom_distance <- 0
top_distance <- max(grid_dm)

# normalizing both training and grid distances
sample_dnorm <- map_dfc(
        sample_dm,
        normalize,
        bottom = bottom_distance,
        top = top_distance
)

grid_dnorm <- map_dfc(
        grid_dm,
        normalize,
        bottom = bottom_distance,
        top = top_distance
)

# this isn't z, it's no
sample_no2norm <- normalize(
        train_coords$CalibratedNO2, 
        top = max(train_coords$CalibratedNO2), 
        bottom = min(train_coords$CalibratedNO2)
)

# and now we can interpolate
# "three ML algorithms I perceive to be the most popular and promising for the task"
# compared with classical ordinary kriging

# fxn to easily look at results
autoplot_int <- function(result){
        
        mapdat <- grid %>%
                left_join(result)   
        
        ggplot() +
                facet_wrap(~model, nrow = 1) +
                geom_sf(data = mapdat,
                        mapping = aes(#x = x, y = y, 
                                fill = z),
                                col = NA) +
                geom_sf(data = test_sf,
                        aes(fill = CalibratedNO2),
                        pch = 22, col = "red") +
                geom_sf(data = sample_sf,
                        aes(fill = CalibratedNO2),
                        pch = 21) +
                geom_sf(data = aqs_sf,
                        aes(fill = epa_no2),
                        pch = 24) +
                scale_fill_distiller(palette = "Spectral") +
                theme_void() + theme(legend.position = "bottom")
        # to do: add in "test" data
        
}

# idw as an actual baseline
fit_idw <- gstat(formula = CalibratedNO2 ~ 1,
                 data = sample_sf,
                 set = list(idp = 2))

# predict across the grid
pred_idw <- predict(fit_idw, grid_sf) %>%
        mutate(z = var1.pred, model = "IDW") %>%
        st_drop_geometry() %>%
        dplyr::select(z, model)
pred_idw$u <- grid_sf$u

# kriging as the baseline

res_krige <- automap::autoKrige(
        formula = CalibratedNO2 ~ 1,
        input_data = as_Spatial(sample_sf),
        new_data = as_Spatial(grid_sf))

pred_krige <- res_krige$krige_output %>%
        as_tibble() %>%
        mutate(z = var1.pred,
               model = "Kriging") %>%
        dplyr::select(z, model)
pred_krige$u <- grid$u

# look at results:        
autoplot_int(rbind(pred_idw, pred_krige))

# neural network: powerful, mostly opaque
# input data: sample point distances and outcomes at sample points
# all normalized
# "train"
nn_input <- cbind(
        z = sample_no2norm,
        sample_dnorm
)
mod_nn <- neuralnet(
        z ~., # col z is dependent var; all others are ind vars
        data = nn_input,
        hidden = c(
                200, 100, 50, 150
        ) # four layers of several arbitrarily chosen neurons
        # what function is getting used in the final transformations
)
# predict for each grid point
pred_nn <- predict(mod_nn, grid_dnorm) %>%
        as_tibble()
# clean tibble of results
pred_nn$z <- denormalize(
        pred_nn$V1, 
        top = max(train_coords$CalibratedNO2), 
        bottom = min(train_coords$CalibratedNO2)
)
pred_nn <- pred_nn %>%
        mutate(model = "Neural Network") %>%
        dplyr::select(z, model)
pred_nn$u <- grid$u

autoplot_int(rbind(pred_idw, pred_nn))        


# SVMs
# "subdivide input data" - draw a convex hull and separating polygons in dims
# require choice of an adequate kernel algorithm
mod_svm <- ksvm(
        z ~ .,
        data = nn_input,
        kernel = "polydot", # this is the decision point
        C = 25 # A parameter to penalize overfitting
)

pred_svm <- predict(mod_svm, grid_dnorm) %>%
        as_tibble()
pred_svm$z <- denormalize(
        pred_svm$V1, 
        top = max(train_coords$CalibratedNO2), 
        bottom = min(train_coords$CalibratedNO2)
)
pred_svm <- pred_svm %>%
        mutate(model = "Support Vector Machine") %>%
        dplyr::select(z, model)
pred_svm$u <- grid$u

autoplot_int(pred_svm)

# random forests
# "a parliament of if-else-clauses and linear models in a trenchcoat"
# params: number of decision trees; no. variables considered each time
mod_rf <- ranger(z~.,
                 data = nn_input, 
                 num.trees = 1000, mtry = 50)

pred_rf <- predict(mod_rf, grid_dnorm)$predictions %>%
        as_tibble()

pred_rf$z <- denormalize(
        pred_rf$value, 
        top = max(train_coords$CalibratedNO2), 
        bottom = min(train_coords$CalibratedNO2)
)

pred_rf <- pred_rf %>%
        mutate(model = "Random Forest") %>%
        dplyr::select(z, model)
pred_rf$u <- grid$u

autoplot_int(pred_rf)

# compare all results:
autoplot_int(rbind(pred_idw, pred_krige,
                   pred_nn, pred_svm,
                   pred_rf))

#### evaluate ------------------------------------------------------------------

final <- rbind(# citywide_mean,
                pred_idw, pred_krige,
               pred_nn, pred_svm,
               pred_rf) %>%
       left_join(st_drop_geometry(test_sf) %>%
                         dplyr::select(#date_day, 
                                 u, y_test = CalibratedNO2)) %>%
        left_join(st_drop_geometry(aqs_sf) %>%
                          dplyr::select(#date_day, 
                                  u, y_val = epa_no2))

eval_interp <- function(y, yhat){
        n <- length(y)
        
        if(n != length(yhat)){
              return(print("Different lengths"))  
        }else{
                # rmse and nrmse
                rmse <- sqrt(sum((y - yhat)^2)/n) 
                nrmse <- rmse / mean(y) # can also be divided by sum(y^2)
                
                #r2
                SSR = sum((y - yhat)^2)
                SST = sum((y - mean(y))^2)
                r2 <- 1 - ((SSR)/(SST))
                
                # rho
                rho <- cor(y, yhat)
                
                return(data.frame(rmse, nrmse, r2, rho, n))
        }
        
}

for(i in unique(final$model)){
        print(i)
        
        eval_interp(final$y_test[!is.na(final$y_test) & final$model == i], 
                    final$z[!is.na(final$y_test) & final$model == i]) %>%
                print()
        
        eval_interp(final$y_val[!is.na(final$y_val) & final$model == i], 
                    final$z[!is.na(final$y_val) & final$model == i]) %>%
                print()
        
        gg <- ggplot() +
                geom_point(data = final[!is.na(final$y_test) & final$model == i,],
                           aes(x = z, y = y_test)) +
                geom_point(data = final[!is.na(final$y_val) & final$model == i,],
                           aes(x = z, y = y_val), col = "red") +
                geom_smooth(data = final[!is.na(final$y_test) & final$model == i,],
                            aes(x = z, y = y_test), 
                            method = "lm")
        
        print(gg)
        
                
}



# RMSE
sqrt(sum((y - yhat)^2)/n)

# NRMSE
rmse / mean(y) # can also be divided by sum(y^2)

# R2
SSR = sum((y - yhat)^2)
SST = sum((y - mean(y))^2)
1 - ((SSR)/(SST))

# MAE
# scatter


#### old, extra code -----------------------------------------------------------

train_day <- train %>%
        # what's the distance
        st_transform(32616) %>%
        # look at one day at a time
        filter(date_day == unique(train$date_day)[1]) %>%
        # keep coords (for universal kriging)
        cbind(st_coordinates(.))

ggplot() + 
        # replace w/ IDW predication
        geom_sf(data = grid) +
        # show points over top
        geom_sf(data = train_day, aes(fill = CalibratedNO2), pch = 21, size = 2) + 
        scale_fill_distiller(palette = "Spectral") +  
        theme_minimal() #+ theme(legend.position = "bottom" ) 


# variogram: function describing relationship between distance and 
# empirical variogram: actual distribution of spatial dependencies in the data
v_emp_ok <- gstat::variogram(
        CalibratedNO2~1,
        as(train_day, "Spatial") # switch from sf to sp
)
plot(v_emp_ok) # hahhhhha oh shit there's no covariance happening?
# alpha and beta should control the correlation
ggplot() + 
        geom_point(data = v_emp_ok, aes(x = dist, y = gamma)) + 
        geom_point(data = v_emp_ok2, aes(x = dist, y = gamma), col = "red", size = 5)

# if/when you are ready:
# (ideally wouldn't autmap...)
v_mod_ok <- automap::autofitVariogram(
        CalibratedNO2 ~ 1,
        as(train_day, "Spatial")
)
plot(v_mod_ok)

ok <- krige(
        CalibratedNO2 ~ 1,
        as(train_day, "Spatial"),
        grid,
        model = v_mod_ok$var_model
)

# simple kriging: mean in the target area is constant and known; local variability is a deviation
# ordinary kriging: mean is constant but unknown
# universal kriging: "not a level plan with some bumps in it but a tilted or curved surface with bumps"

