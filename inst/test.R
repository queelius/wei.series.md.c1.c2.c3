


library(wei.series.md.c1.c2.c3)
library(dplyr)
library(md.tools)


df <- generate_guo_weibull_table_2_data(30)
df <- md_decorate_component_cause_k(df)
res <- wei.series.md.c1.c2.c3::empirical_candidates_k(df)
df
res
