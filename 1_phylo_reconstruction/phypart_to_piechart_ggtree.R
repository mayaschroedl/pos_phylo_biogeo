


# Working directory -------------------------------------------------------
wd = file.path( getwd(), "1_phylo_reconstruction", "5_phypartstopiecharts", "2")



# Phypart input -----------------------------------------------------------
phypart_out_file = file.path(wd, ".hist")

# you cannot directly open the file as a dataframe, because of different numbers of columns per line; one has to do a little tweak
max_col = max(count.fields(phypart_out_file, sep=",")) #maximum number of columns
phypart_out = read.table(file.path(wd, ".hist"), h=F, sep=",", col.names = paste0("V",seq_len(max_col)), fill = TRUE) # not we can read the file as dataframe


phypart_out$node = phypart_out$V1
phypart_out$V1
