#todo: test and describe

GroupSamples <- function(samples.dt, groups.dt)
{
  setkey(samples.dt, sample_id_1)
  setkey(groups.dt, sample_id)
  group_ids <- as.vector(unique(groups.dt[, group_id]))
  
  joined.dt <- samples.dt[groups.dt, nomatch=0]
  joined.dt[, group_ids[1]:= ifelse(group_id==group_ids[1], sample_id_1, sample_id_2)]
  joined.dt[, group_ids[2]:= ifelse(group_id==group_ids[2], sample_id_1, sample_id_2)]
  result.dt <- joined.dt[, c(group_ids[1], group_ids[2]), with=F]
  return(result.dt)
}


make.short.name <- function(vector.names) {
  #get specie from new names if OTU
  vector.names.new <- str_replace(vector.names, ".*NN=", "")
  vector.names.new <- str_replace(vector.names.new, "\\|D=.*", "")
  #for new names - clear unclassified
  vector.names.new <- str_replace_all(vector.names.new, ";unclassified", "")
  #for old names - clear unclassified
  vector.names.new <- str_replace_all(vector.names.new, " .__$| .__;", "")
  vector.names.new <- str_replace_all(vector.names.new, "__", ":")
  vector.names.new <- str_replace_all(vector.names.new, ";$", "")
  # get for all two last entities
  vector.names.new <- gsub('.*;\\s*([^;]+);\\s*([^;]+)$', '\\1 \\2', vector.names.new, perl=T)
  vector.names.new <- strsplit(vector.names.new, "[; _]+")
  vector.names.new <- sapply(vector.names.new, function(x) {
    ifelse(length(unique(x))<=2, paste(unique(x), collapse=' '), paste(unique(x)[1], unique(x)[2], collapse=' '))})
  return(vector.names.new)
}

aggregate.by.name <- function(data) {
  nms <- colnames(data)
  data1 <- sapply(unique(nms), function(i) {
    if(is.matrix(data[, nms==i])) {rowSums(data[, nms==i])
      } else {data[, nms==i]}
    })
  return(data1)
}
