flog.info("TOTAL START")

################################################
# NeededSampleSize
################################################
source("functions_from_Vera/functions_power.R")

#loading data
genus.bef <- genus[rownames(genus) %in% data.samples$after,]
genus.aft <- genus[rownames(genus) %in% data.samples$before,]

#собираем малопредставленные виды в один столбец, умножаем сначала данные на константу, 
#чтобы все ненулевые элементы матрицы были больше 1
#ToDo: сделать сранение всех групп из набора между собой
filt.before <- Data.filter(genus.bef*6000, "data", 1, 30)
filt.after <- Data.filter(genus.aft*6000, "data", 1, 30)

data <- merge.features(filt.before,filt.after, 0)
group.before <- data$dataCase
group.after <- data$dataCntrl

power.curve(group.before, group.after, 1 ,2000)
flog.info("END")
