library(haven)

fns <- grep("sas7bdat", list.files("data/uat_dataset_extract"), value = TRUE)
names <-  sub("\\.sas7bdat", "", fns)
paths <- file.path("data/uat_dataset_extract", fns)

for (i in seq_len(length(names))) {
    assign(names[i], read_sas(paths[i]), envir = .GlobalEnv)
}
