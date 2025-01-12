library(readxl)
library(openxlsx)
library(dplyr)
library(tidyverse)
library(purrr)


setwd("C:/Users/yuqian.wang/NTU_Sherry/8Mpox/submit")
path <- "C:/Users/yuqian.wang/NTU_Sherry/8Mpox/data/ID.xlsx"

#setting
DL <- 3 # Detection limit

# Get the sheet names
sheet_names <- excel_sheets(path)

# Read each sheet into a named list
original <- setNames(lapply(sheet_names, function(sheet) read_excel(path, sheet = sheet)), sheet_names)


# Count rows for each ID
calculate_sample <- function(original){
  
  id_counts <- original %>%
    group_by(ID) %>% 
    summarise(RowCount = n(), .groups = "drop")
  # Count how many IDs have 3 data points, and how many have 2
  count_VL <- id_counts %>%
    count(RowCount)
  
  return(count_VL)
}

count_VL<-map(original,calculate_sample)

combine_count <- map_df(count_VL, ~as.data.frame(.x),.id = "site")

combine_VL <- map_df(original[c('Rectum','Saliva','Oropharynx')], ~as.data.frame(.x),.id = "site") %>%
  mutate(VL = round(as.numeric(ifelse(`log（copies/ML）` == "Negative", 3, (`log（copies/ML）`))),2),
         censor = ifelse(`log（copies/ML）` == "Negative", 1, 0)) %>%
  group_by(ID,site) %>%
  filter(n() >= 1) %>%  # Keep only IDs that appear more than twice
  ungroup() %>%
  unite(ID_site, c("ID", "site"),remove = FALSE)

write.csv(combine_VL, file = 'data/combine_VL_site_1point.csv',row.names=FALSE)


# Split the data frame by a specific column, say 'category'
split_VL <- split(combine_VL, combine_VL$site)

output_path <- "data/" 

# Save each split data frame to a separate Excel file
for (name in names(split_VL)) {
  # Create a new workbook for each split data frame
  wb <- createWorkbook()
  # Add the split data to the workbook
  addWorksheet(wb, "Sheet1")
  writeData(wb, sheet = "Sheet1", split_VL[[name]])
  
  # Save the workbook to a different Excel file, named by the category
  file_name <- paste0(output_path,name, "_VL_3", ".xlsx")
  saveWorkbook(wb, file_name, overwrite = TRUE)
}



#Combine pre-symptomatic and symptomatic data#####
pre_sym <- read_csv("data/Mpox data.csv")

#unique value
unique(pre_sym$`Sample Location`)
unique(combine_VL$site)

sub_presym <- pre_sym %>% 
  filter(`Sample Location` %in% c('Anorectal Swab','Oropharyngeal Swab','Saliva')) %>%
  mutate(site = case_when(`Sample Location` == "Anorectal Swab" ~ "Rectum",
                          `Sample Location` == "Oropharyngeal Swab" ~ "Oropharynx",
                          `Sample Location` == "Saliva" ~ "Saliva"),
         `Days post symptoms onset` = Day - systematic_symptom_onset,
         raw_VL = (`Ct Value`-41.388)/-3.611,
         VL = ifelse(raw_VL <= 3, 3, raw_VL),
         censor = ifelse(VL <= 3, 1, 0)) %>%
  unite(ID_site, c("Patient", "site"),remove = FALSE)

ggplot(sub_presym) +
  geom_point(aes(x=`Days post symptoms onset`,y=VL, colour=site)) +
  geom_line(aes(x=`Days post symptoms onset`,y=VL, colour=site)) +
  facet_wrap(vars(Patient), ncol=10, nrow=8)

#combine before and after symptom onset
common_columns <- intersect(names(combine_VL), names(sub_presym))
df1_common <- combine_VL[, common_columns, drop = FALSE]
df2_common <- sub_presym[, common_columns, drop = FALSE]

combine_sym <- rbind(df1_common, df2_common)
write.csv(combine_sym, file = 'data/combine_VL_pre_and_sym_1point_limit3.csv',row.names=FALSE)
  
  