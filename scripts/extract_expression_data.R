library(tidyr)
library(dplyr)
df_ = data.frame(id = c(), gene = c(), expr = c(), age = c(), age_cat = c(), sex = c(), stim = c())
for (stim in c("null","ecoli","bcg","saureus","seb","candida","flu")) {
  print (stim)
  prep = as.data.frame(t(data[[stim]]$expression))
  prep$id = as.integer(rownames(prep))
  prep = dplyr::left_join(prep, age_sex[,c("id","sex","AGE.V0","AGE.Cat")], by = "id")
  prep = tidyr::gather(prep, "Gene", "Expression", ABCB1:TBP)
  colnames(prep) = c("id", "sex", "age", "age_cat", "gene", "expr")
  prep$stim = stim
  prep = prep[,c("id","gene", "expr", "age", "age_cat", "sex", "stim")]
  batch_info = data[[stim]]$info[,c("SUBJID","batch.TruCult","batch.NS")]
  prep = dplyr::left_join(prep, batch_info, by = c("id"="SUBJID"))
  colnames(prep) = c("id","gene", "expr", "age", "age_cat", "sex", "stim", "truculture_batch", "ns_batch")
  df_ = dplyr::bind_rows(df_, prep)  
}
save(df_, file="./data/gene_response_age.corrected.Rdata")
head(df_)
