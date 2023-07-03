library(dplyr)
library(readr)

rootdir <- '~/project/CATD-ReliabilityAnalysis/'
datadir <- paste0(rootdir, 'data/')

#Set this to home directory
rootdir <- '~/project/CATD-ReliabilityAnalysis/'
#Subject-level information
all_subs_unfiltered <- read_csv(paste0(rootdir,'references/rest_df.csv'), progress = TRUE)

#tanner pubertal data
tanner <- read_csv(paste0(rootdir, 'references/tanner.csv'))

# #Demographic information
dm <- read_csv(paste0(rootdir,"references/rest_paths_no_dates_race.csv"))


all_subs <- merge(all_subs_unfiltered,
                  dm[,c("subject","session","p_demo_screen_background_race")], 
                  by=c("subject","session"),all.x=TRUE)

all_subs <- merge(all_subs, tanner, by.x = c("subject","session"), by.y = c("SDAN", "Task_Visit_Type"), all.x = T)

all_subs$old_sess <- all_subs$session
all_subs$session[grepl("v1",all_subs$session)] <- "v1"
all_subs$session[grepl("v4",all_subs$session)] <- "v4"
all_subs$session[grepl("v2",all_subs$session)] <- "v2"
all_subs$session[grepl("v3",all_subs$session)] <- "v3"

write_csv(all_subs, paste0(rootdir, 'references/all_subjects_df.csv'))
