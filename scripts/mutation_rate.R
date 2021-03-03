Args <- commandArgs()

reformat = read.delim(Args[6],header = F,stringsAsFactors = F)
reformat <- reformat[!duplicated(paste(reformat$V1,reformat$V2)),]
reformat$V2 <- gsub("reversed\\|","",reformat$V2)

mutation_rate <- function(x)
{
  n = nrow(x)
  b = x$V9[n] - x$V8[1] + 1
  t = sum(x$V6) + sum(x$V7) + b - sum(x$V5)
  return(t/b)
}

mr <- sapply(split(reformat,reformat$V2), mutation_rate)
write.table(mr, Args[7], sep = '\t', row.names = T, col.names = F, quote = F)


