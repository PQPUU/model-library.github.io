library(dplyr)

df <- read.table("simtab.dat", header = TRUE) %>% #reading in the simulated data, one record per ID
  select(-RTTE, -SURX, -ICOUNT, -ITER, -RAND) #select only needed variables

df <- rbind(df %>% mutate(TIME = 0, #initiation record with time zero
                          TYPE = 1, #internal filtering variable
                          EVID = 3, #event ID
                          DV = 0),  #dependent variable representing no event at time 0
                 
                 df %>% mutate(TIME = ifelse(DV==1, ceiling(TIME), TIME), #rounding of time to get days, actual time event/censoring
                                    TYPE = 2, #internal filtering variable
                                    EVID = 0), #event ID
                 
                 df %>% mutate(TIME = max(ENDTIME), #maximum time for simulations
                                    TYPE = 3, #internal filtering variable
                                    EVID = 0,#event ID
                                    DV = 0)) %>% #dependent variable representing no event, to be overwritten in simulations
  
  arrange(ID, TYPE) %>% #arrangint the dataset based on ID and time
  rename(FTIME = ENDTIME) %>% #renaming ENDTIME to FTIME
  select(ID, DV, TIME, EVID, TYPE, FTIME, GROUP) #selecting desired columns

#creating a dataset
write.csv(df, "Simulated_dataset_TTE_hospitalization_BCG_only.csv", quote = F, row.names = F)