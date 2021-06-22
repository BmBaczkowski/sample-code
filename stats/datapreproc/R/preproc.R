# --- 
library("tidyverse")
# ---

# --- specify data directory
dat_dir       <- file.path("..", "data", "jatos_resultfiles")
# --- read json file with the key of ids
json_key      <- jsonlite::read_json(file.path(dat_dir, "key.json"))
# --- load csv with prolific data 
prolific      <- read_csv(file.path("..", "data", "prolific_export_day01.csv"))
# --- create a list of subject ids
dir_list      <- dir(dat_dir, pattern = "[0-9]")
# --- create a list of csv files with memory task
memory_list   <- lapply(dir_list, function(n) {
                              fname <- file.path(dat_dir, n, "memory.csv")
                              subid <- n
                              return(c(subid, fname))
                              })
# --- create a list of csv files with choice task
choice_list   <- lapply(dir_list, function(n) {
                              fname <- file.path(dat_dir, n, "choice.csv")
                              subid <- n
                              return(c(subid, fname))
                              })
# --- create data frame for memory task
df_memory     <- lapply(memory_list, function(x) {
                        if (file.exists(x[2])) {
                          subid <- names(which(json_key == x[1]))
                          age <- prolific$age[prolific$participant_id == subid]
                          sex <- prolific$Sex[prolific$participant_id == subid]
                          df <- read_csv(x[2])
                          df <- add_column(df, id = rep(x[1], nrow(df)), 
                                               age = rep(age, nrow(df)),
                                               sex = rep(sex, nrow(df)),
                                               .before = "trialid")
                          return(df)
                        }
                      })
df_memory       <- bind_rows(df_memory)
# --- save data frame
saveRDS(df_memory, file.path("..", "data", "df_memory.Rds"))
write_csv(df_memory, path = file.path("..", "data", "df_memory.csv"))
# --- create data frame for choice task
df_choice       <- lapply(choice_list, function(x) {
                        if (file.exists(x[2])) {
                          subid <- names(which(json_key == x[1]))
                          age <- prolific$age[prolific$participant_id == subid]
                          sex <- prolific$Sex[prolific$participant_id == subid]
                          df <- read_csv(x[2])
                          df <- add_column(df, id = rep(x[1], nrow(df)), 
                                           age = rep(age, nrow(df)),
                                           sex = rep(sex, nrow(df)),
                                           .before = "trialid")
                          return(df)
                        }
                      })
df_choice         <- bind_rows(df_choice)
# --- save data frame
saveRDS(df_choice, file.path("..", "data", "df_choice.Rds"))
write_csv(df_choice, path = file.path("..", "data", "df_choice.csv"))

