#########################################################################
# Compare fold changes of selected taxa between the groups.
# The taxa that get MAX_RATIO are those that "appear" or "disappear".
#########################################################################

# some good threshold presets
MIN_CHANGE_RATIO <- 1.5 # minimum fold change (e.g. increase OR decrease )
MIN_MAX_PERC <- 1 # minimum threshold for maximum of taxon rel. abundance across 2 time points

# let's test for one subject
tag_before <- data.samples$before[10]
tag_after <- data.samples$after[10]
#tag_before <- data.samples$before[1]
#tag_after <- data.samples$after[1]
#tag_before <- data.samples$before[2]
#tag_after <- data.samples$after[2]

# families
ff <- family[c(data.samples$before, data.samples$after),]
rs <- get_changing_bugs_for_a_subject(ff, tag_before, tag_after, MIN_MAX_PERC, MIN_CHANGE_RATIO)
rs
rs$taxa_decreased
rs$taxa_increased

# genera
ff <- genus[c(data.samples$before, data.samples$after),]
rs <- get_changing_bugs_for_a_subject(ff, tag_before, tag_after, MIN_MAX_PERC, MIN_CHANGE_RATIO)
rs

# species
ff <- species[c(data.samples$before, data.samples$after),]
rs <- get_changing_bugs_for_a_subject(ff, tag_before, tag_after, MIN_MAX_PERC, MIN_CHANGE_RATIO)
rs
