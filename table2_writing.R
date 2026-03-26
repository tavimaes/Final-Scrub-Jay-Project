# ----- combining stats for table 2 ------

binary_data <- read_csv("binary_behav_stats.csv")

call_data <- read.csv("call_count_stats.csv")

table2 <- bind_rows(binary_data, call_data)


write_csv(table2, "table2.csv")
