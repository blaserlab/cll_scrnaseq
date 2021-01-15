source("00_packages_functions.R")

ind_qc_res <- pmap(.l = list(cds = cds_list, 
                             genome = rep("human", times = length(cds_list))),
                   .f = function(cds, genome) {
                     return_list <- qc_func(cds = cds, genome = genome)
                     return(return_list)
                   })
names(ind_qc_res) <- names(gex_pipestance_list)

# inspect the plots
ind_qc_res[[4]][[1]]
# ind_qc_res[[2]][[3]]
# ind_qc_res[[3]][[3]]
# ind_qc_res[[4]][[3]]
