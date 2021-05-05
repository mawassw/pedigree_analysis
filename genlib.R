###########################@ GENLIB functions @############################
library("GENLIB")
?`GenlibR-package`

saguenay <- as.data.frame(read_excel("C:/Users/walid/OneDrive/Bureau/Work/Pedigree/Saguenay_ind.xlsx"))

saguenay[,6] <- NULL

saguenay$datn <- as.numeric(format(as.Date(saguenay$datn, format = "%Y-%m-%d"),"%Y")) #extract only birth year

saguenay <- saguenay[saguenay$sexe != 9,] #remove individuals with unknown sex

####setting up genlib object
saguenay_genlib <- cbind.data.frame("ind"=saguenay$ind,"father"=saguenay$pere,"mother"=saguenay$mere,"sex"=saguenay$sexe)

sag_gen <- gen.genealogy(saguenay_genlib, autoComplete = TRUE)

####getting number of children
ped_sag$n_child <- rep(0,nrow(ped_sag))
ind <- ped_sag$ind

for (i in ind) {
  ped_sag[ped_sag$ind == i,]$n_child <- gen.nochildren(sag_gen, individuals = i)
}