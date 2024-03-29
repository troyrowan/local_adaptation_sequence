

get_animal_table <- function(date_string){
#animal_table <- Hmisc::mdb.get(str_c("/Volumes/UMAG_samplesDB/Samples_", date, ".accdb"), tables = "Animal", allow = c("_"))
df <- Hmisc::mdb.get(glue::glue("~/Desktop/Samples_{date_string}.accdb"), tables = "Animal", allow = c("_"))
  
  #animal_table <- Hmisc::mdb.get(glue::glue("/Users/harlyjanedurbin/Box Sync/HairShedding/ReportedData/Samples_190907.accdb"), tables = "Animal", allow = c("_"))
  
df <- doctoR::clear_labels(df) %>% 
  dplyr::mutate_if(is.factor, as.character) 

#Tried to do this in a pipable fashion using the naniar package but the data frame is too large
df[df == ""] <- NA

#for now, append AMGV to animal table Gelbvieh

df <-
  df %>% 
  # Remove white space
  dplyr::mutate_if(is.character, funs(stringr::str_squish(.))) %>% 
  dplyr::mutate(breed_assoc =
                  dplyr::case_when(
                    breed_assoc %in% c("American Angus Association ", "American Anugs Association") ~ "American Angus Association",
                    breed_assoc %in% c("American Shorthorn Association"," American Shorthon Association") ~ "American Shorthorn Association",
                    breed_assoc %in% c("American Maine Anjou Association", "Maine Anjou Association of America") ~ "American Maine-Anjou Association",
                    TRUE ~ as.character(breed_assoc)
                  ),
                Comment = dplyr::case_when(
                  Comment %in% c("Local Adaptation Project ", "Local Adptation Project", "Local Adaptaion Project", "Local Adapatation Project") ~ "Local Adaptation Project", 
                  TRUE ~ as.character(Comment)
                ),
                Reg = dplyr::case_when(
                  BC == "GEL" & breed_assoc == "American Gelbvieh Association" & breed_assoc_country == "USA" & !stringr::str_detect(Reg, "GV") ~ stringr::str_c("AMGV", Reg), 
                  TRUE ~ as.character(Reg)
                )
  ) %>% 
  dplyr::mutate(DOB = stringr::str_remove(DOB, " 00:00:00"),
         DOB = lubridate::mdy(DOB))

return(df)
}