prooject.dir <- '/projects/16s_study'
setwd(prooject.dir)
library(data.table)
library(stringr)

data <- fread('data_2016_04_28/meta_tables/biotabarcodes.csv')
data$posix <- as.numeric(strptime(data$created, "%d.%m.%Y %H:%M:%S"))
data[,paired:=ifelse(.N==2, TRUE, FALSE),by=user_id]
before <- data[paired==TRUE, .SD[which.min(posix)], by=user_id]
after <- data[paired==TRUE, .SD[which.max(posix)],by=user_id]
before$barcode <- str_replace(before$barcode, '-', '.')
after$barcode <- str_replace(after$barcode, '-', '.')

#PAIRED SAMPLES
data.samples <- data.table(cbind(before$barcode, after$barcode))
setnames(data.samples, c("sample_id_1", "sample_id_2"))

#GROUPS
data.groups.1 <- data.table(sample_id = data.samples$sample_id_1, group_id= 'BEFORE_TREAT', quality_status='good')
data.groups.2 <- data.table(sample_id = data.samples$sample_id_2, group_id= 'AFTER_TREAT', quality_status='good')
data.groups <- rbindlist(list(data.groups.1, data.groups.2))

#RUNS
#runs <- data.table(sample = data.samples$sample_id_1 + data.samples$sample_id_2, run_id = 1)
#RUNS_DESCRIPTION
#runs <- data.table(sample = data.samples$sample_id_1 + data.samples$sample_id_2, run_id = 1)

write.table(data.samples, file='data_2016_04_28/meta_tables/16s_pairedsamples.csv', sep=",", quote = F, row.names = F)
write.table(data.groups, file='data_2016_04_28/meta_tables/16s_samples.csv', sep=",", quote = F, row.names = F)
