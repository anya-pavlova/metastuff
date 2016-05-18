write.facts <- function(tax.matrix, dir.sample.id, fact.table, bact.catalog, cut.names=F, append=F) {
  dict <- data.table('local_name'=colnames(tax.matrix))
  dict$name.merge <- dict$local_name
  dict$name.merge <- str_replace_all(dict$name.merge, '\\w:', '')
  if (cut.names) {
    dict[, name.merge := str_split(name.merge, ' ')]
    dict$name.merge  <- sapply(dict$name.merge , function(x) ifelse(length(x)>1, x[2], x[1]))
  }
  dict <- merge(dict, bact.catalog, by.x='name.merge', by.y = 'name', all.x=T)
  
  n <- round(nrow(dict[!is.na(id)])*100/nrow(dict), 1)
  message(n, ' % of bacteria from DataBase were found in NCBI.')
  
  dict <- dict[!is.na(id)]
  local.facts <- merge(fact.table, dict, by='id')
  local.facts <- local.facts[, c('from', 'fact', 'to', 'id', 'dis_id', 'local_name', 'N'), with=F]
  levels=apply(tax.matrix[,colnames(tax.matrix) %in% local.facts$local_name], 2, calculate.level)
  rownames(levels) <- rownames(tax.matrix)
  
  dir.sample.id[, {
    dt.list <- lapply(c('low', 'high'), function(x){ 
      facts.before <- get.facts(before, x, levels, local.facts, 'BEFORE')
      facts.after <- get.facts(after, x, levels, local.facts, 'AFTER')
      return(rbindlist(list(facts.before, facts.after)))
    })
    dt <- rbindlist(dt.list)
    dt <- setnames(dt, c('N'), c('score'))
    dt <- setorderv(dt, c('group', 'score'), c(-1, -1))
    dt <- setcolorder(dt, c('group', 'abundance', 'score', 'from', 'fact', 'to', 'id', 'dis_id', 'local_name', 'sample_id'))
    write.table(dt, file.path(dir, paste0(before, '_', after, 'disease_facts.csv')),
                sep='\t', col.names = T, row.names = F, quote = F, append=append)
  }, by=1:nrow(data.samples)]
}

get.facts <- function(id, level.name, levels.data, facts.data, group.name) {
  i.bact <- colnames(levels.data)[levels.data[id,]==level.name]
  i.facts <- facts.data[local_name %in% i.bact]
  i.facts$group <- group.name
  i.facts$sample_id <- id
  i.facts$abundance <- level.name
  return(i.facts)
}
