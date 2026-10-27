library(scTenifoldKnk)
for(n in c('scTenifoldKnk','knockoutDirection')) {cat('\nFUNCTION ',n,'\n');print(get(n,asNamespace('scTenifoldKnk')))}