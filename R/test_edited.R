# Title     : TODO
# Objective : TODO
# Created by: sarahleinicke
# Created on: 3/14/18

# Reference cell clustering script
source('celda_functions.R')
source('celda.R')
source('celda_C.R')
source('celda_G.R')
source('celda_CG.R')
source('split_clusters.R')
source('s3_generics.R')
library(doParallel)
#library(parallel)
#source("https://bioconductor.org/biocLite.R")
#biocLite("BiocParallel")
#library("BiocParallel")
#library(celda)
#library(hashmap)

main = function(){
    #time_celda();
    #time_celda_CG_simCells();
    time_celda_CG();
    #time_celda_C();
    #time_celda_G();
}


###############
#    CELDA    #
###############

time_celda = function(){

    print('CELDA - celdaCGsim.rda data')
    load("../tests/celdaCGsim.rda")
    print('Dimensions - celdaCGsim.rda')
    print(dim(celdaCG.sim$counts))

    print('EDITED')
    print(system.time(celda(counts = celdaCG.sim$counts,
        model = "celda_CG", nchains = 1, K = celdaCG.sim$K, L = celdaCG.sim$L,
        max.iter = 15, random.state.order = FALSE)))

}



#################
#    CELDA CG   #
#################

time_celda_CG_simCells = function(){
    # Default simulateCell ranges:
    #    C.Range=c(50,100)   - cell counts
    #    N.Range=c(500,5000) - transcript counts
    #    G = 1000            - num genes

    print('CELDA_CG - simulateCells.celda_CG data')
    #load("large_simc_data.rda");
    sim_counts_CG = simulateCells.celda_CG(
        "celda_CG", S=10, C.Range=c(50, 100), N.Range=c(500,5000), G=1000, K=5, L=9)

    print('Dimensions - simulateCells_celda_CG')
    print(dim(sim_counts_CG$counts))

    print('EDITED')
    print(system.time(celda_CG(counts = sim_counts_CG$counts,
        K = 5, L=9, max.iter = 15,
        random.state.order=FALSE,)))
}


time_celda_CG = function(){
    print('CELDA_CG - celdaCGsim.rda data')
    load("../tests/celdaCGsim.rda")

    print('Dimensions - celdaCGsim.rda')
    print(dim(celdaCG.sim$counts))

    print('EDITED')
    print(system.time(celda_CG(counts = celdaCG.sim$counts,
        K = celdaCG.sim$K, L=celdaCG.sim$L, max.iter = 15,
        random.state.order = FALSE)))
}



#################
#    CELDA C    #
#################

time_celda_C = function(){
    print('CELDA_C')
    sim_counts_C = simulateCells.celda_C()
    print('Dimensions - simulateCells.celda_C()')
    print(dim(sim_counts_C$counts))

    print('EDITED - with parApply')
    print(system.time(celda_C(counts = sim_counts_C$counts,
        sample.label=NULL, K = 5, alpha=1, beta=1,
        stop.iter = 10, max.iter=15, split.on.iter=10, split.on.last=TRUE,
        random.state.order=FALSE, count.checksum=NULL, seed=12345,
        z.init = NULL, process.counts=TRUE, logfile=NULL)))

}

#################
#    CELDA G    #
#################


time_celda_G = function(){
    print('CELDA_G')
    #sim_counts_G = simulateCells.celda_G()
    #save(sim_counts_G, file='sim_counts_G.rda')

    load('sim_counts_G.rda')
    print('Dimensions - simulateCells.celda_G()')
    print(dim(sim_counts_G$counts))

    print('EDITED')
    print(system.time(celda_G(counts = sim_counts_G$counts,
        L = 5, beta=1, gamma = 1,
        stop.iter = 10, max.iter= 200, split.on.iter=10, split.on.last=TRUE,
        random.state.order=FALSE, count.checksum=NULL, seed=12345,
        y.init = NULL, process.counts=TRUE, logfile=NULL)))
}



main();





