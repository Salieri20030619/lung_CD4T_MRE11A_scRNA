cat('DR knk exported:', 'dRegulation' %in% getNamespaceExports('scTenifoldKnk'),'\n')
cat('DR net exported:', 'dRegulation' %in% getNamespaceExports('scTenifoldNet'),'\n')
cat('MA net exported:', 'manifoldAlignment' %in% getNamespaceExports('scTenifoldNet'),'\n')
print(scTenifoldNet::makeNetworks)