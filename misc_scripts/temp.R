x = lapply(1:18, function(x) list())
y = read.table("reclass.age_covar.txt", stringsAsFactors = F, sep = "\t", header = T)
use_names = unique(y[,2])

i <- 1
for(f in list.files(pattern="reclass")){
  print(i)

  y = read.table(f, stringsAsFactors = F, sep = "\t", header = T)
  for(z in 1:length(use_names)){
    x[[z]][[i]] = y[y[,2] == use_names[z],]
  }
  if(i == 5){
    for(z in 1:length(use_names)){
      x[[z]] <- do.call("rbind", x[[z]])
    }
  }
  i <- i + 1
}




