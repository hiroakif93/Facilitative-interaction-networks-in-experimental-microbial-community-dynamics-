####################
# -- infulence code
# -- Created by K.Suzuki

## functions
nm <- function(x, z){
    x.z = cbind(abs(x),z)
    x.sum = tapply(x.z[,1], x.z[,2], sum)[x.z[,2]]
    xm = x.z[,1] / as.numeric(x.sum)
    xm}

influence <- function(mat, mode){
    l=nrow(mat)

    if (mode=='all') {
        matz <- abs(mat)
    } else if (mode=='positive') {
        mat[mat<0] <- 0
        matz <- mat
    } else if (mode=='negative') {
        mat[mat>0] <- 0
        matz <- abs(mat)
    } else {
        print("mode: Only can use \"all\", \"positive\" or \"negative\"")
    }

    I <- diag(l)
    matz=matz+0.001*I
    mm=matz
    mm[mm>0] <- 1
    g=graph_from_adjacency_matrix(mm)
    mem=components(g, mode = "weak")[['membership']]

    x.z = cbind(rep(1,length(mem)), mem)
    mnum = as.numeric(tapply(x.z[,1], x.z[,2], sum)[x.z[,2]])
    mnum[mnum==1] <- 0
    mnum[mnum!=0] <- 1

    vi <- matrix(nm(rnorm(l), mem))

    ki <- colSums(matz) ##
    ma <- matz/ki-I

    ndv = Inf

    t=0
    nmon <-c()
    while (ndv>0.01 && t <5000) {
    t = t + 1
    dv = ma %*% vi
    ndv = norm(dv*mnum)
    nmon = c(nmon, ndv)
    vi = nm(vi + 0.01 * dv, mem)
    }

    if (ndv>0.01) {print("Convergence failure: vi may incorrect")}
    
    return(matrix(c(vi, mem), l))
}

##################################################