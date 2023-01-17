unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()

# Run ALRA
ALRA.EnImpute2 = function(count, k = 0, q = 10){

  count.t = t(count)
  # library and log normalization
  count_norm = normalize_data(count.t)
  # Impute using the function alra
  count.ALRA = alra(count_norm,k=k, q=q)[[2]]
  count.ALRA = t(count.ALRA)

  count.ALRA = exp(count.ALRA)- 1

  row.names(count.ALRA) = row.names(count)
  colnames(count.ALRA) = colnames(count)
  count.ALRA = as.matrix(count.ALRA)
  count.ALRA
}

# Run DCA
DCA.EnImpute2= function(count, normtype = "zheng", type = "zinb-conddisp",
                        l2 = 0, l1 = 0, l2enc = 0, l1enc = 0, ridge = 0,
                        gradclip = 5, activation = "relu", hiddensize = "64,32,64",
                        hyper = FALSE, hypern = 1000){

  count = round(count)
  dir.create("./DCA_result")
  utils::write.csv(count,"./DCA_result/count.csv")

  comand = paste("dca ./DCA_result/count.csv ./DCA_result/result",
                 "--normtype", eval(normtype),
                 "--type", eval(type),
                 "--l2", eval(l2),
                 "--l1", eval(l1),
                 "--l2enc", eval(l2enc),
                 "--l1enc", eval(l1enc),
                 "--ridge", eval(ridge),
                 "--gradclip", eval(gradclip),
                 "--activation", eval(activation),
                 "--hiddensize", eval(hiddensize))
  if (hyper==TRUE){
    comand =  paste(comand, "--hyper", "--hypern", eval(hypern))
  }

  # Impute using DCA
  system(eval(comand))


  count.DCA = utils::read.csv("./DCA_result/result/mean.tsv", sep="\t",header = TRUE, row.names = 1)
  unlink("./DCA_result", recursive=TRUE)

  count.DCA = as.matrix(count.DCA)
  count.DCA
}

# Run DrImpute
DrImpute.EnImpute2 = function(count, ks = 10:15, dists = c("spearman", "pearson"), method = "mean",
                              cls = NULL){
  # Preprocess gene expression matrix
  count = DrImpute::preprocessSC(count, min.expressed.gene = 0, min.expressed.cell = 0, max.expressed.ratio = 1,
                                 normalize.by.size.effect = TRUE)
  # log transformation
  logcount = log(count+1)
  # Impute using the function DrImpute
  count_DrImpute = DrImpute::DrImpute(logcount, ks = ks, dists = dists, method = method,
                                      cls = cls)

  count_DrImpute = exp(count_DrImpute)-1

  row.names(count_DrImpute) = row.names(count)
  colnames(count_DrImpute) = colnames(count)
  count_DrImpute = as.matrix(count_DrImpute)
  count_DrImpute
}

#Run knnsmooth
knnsmooth.EnImpute2=function(count,k=10){
  count_knnsmooth <- knn_smoothing(count,k=k)

  row.names(count_knnsmooth) = row.names(count)
  colnames(count_knnsmooth) = colnames(count)
  count_knnsmooth = as.matrix(count_knnsmooth)
  count_knnsmooth
}

# Run MAGIC
MAGIC.EnImpute2 = function(count, k = 10, alpha = 15, t = "auto", npca = 20,
                           t.max = 20, knn.dist.method = "euclidean", n.jobs = 1){
  count.t = t(count)
  # Library size normalization
  count.normalized = Rmagic::library.size.normalize(count.t)
  # Log normalization
  count.log = log(count.normalized + 1)
  # Impute using the function magic
  count_MAGIC = Rmagic::magic(count.log, k = k, alpha = alpha, t = t, npca = npca,
                              t.max = t.max, knn.dist.method = knn.dist.method, n.jobs = n.jobs)

  count_MAGIC = t(as.matrix(count_MAGIC))
  count_MAGIC[count_MAGIC<0]=0
  count_MAGIC = exp(count_MAGIC)-1

  row.names(count_MAGIC) = row.names(count)
  colnames(count_MAGIC) = colnames(count)
  count_MAGIC = as.matrix(count_MAGIC)
  count_MAGIC
}

# Run SAVER
SAVER.EnImpute2 = function(count, do.fast = TRUE, ncores = 1, size.factor = NULL,
                           npred = NULL, null.model = FALSE, mu = NULL){
  # Impute using the function saver
  count_SAVER= SAVER::saver(count, do.fast = do.fast, ncores = ncores, size.factor = size.factor,
                            npred = npred, null.model = null.model, mu=mu)$estimate

  row.names(count_SAVER) = row.names(count)
  colnames(count_SAVER) = colnames(count)
  count_SAVER = as.matrix(count_SAVER)
  count_SAVER
}

# Run scImpute
scImpute.EnImpute2 = function(count, drop_thre = 0.5, Kcluster = 10, labeled = FALSE,
                              labels = NULL, genelen = NULL, ncores = 1){

  dir.create("./scImpute_result")
  saveRDS(count, file ="./scImpute_result/count.rds")
  # Run scImpute using the function scimpute
  out = scImpute::scimpute(count_path = "./scImpute_result/count.rds",
                           infile = "rds",
                           outfile = "rds",
                           type = "count",
                           out_dir = "./scImpute_result/",
                           labeled = labeled,
                           drop_thre = drop_thre,
                           Kcluster = Kcluster,
                           genelen = genelen,
                           ncores = ncores)

  count_scImpute = readRDS("./scImpute_result/scimpute_count.rds")
  unlink("./scImpute_result", recursive=TRUE)

  row.names(count_scImpute) = row.names(count)
  colnames(count_scImpute) = colnames(count)

  count_scImpute = as.matrix(count_scImpute)
  count_scImpute
}

#Run scNPF
scNPF.EnImpute2=function(count,network="context",gamma = 0.5, qt.gene = 0.4, qt.cell = 0.5,
                         nThreads = 1){
  unregister_dopar()
  library('WGCNA')#description里面要写上会自动import进来
  count_scNPF<- scNPF::scNPF.pro(x=count, network=network,gamma = gamma, qt.gene = qt.gene, qt.cell = qt.cell,
                                 nThreads = nThreads)

  row.names(count_scNPF) = row.names(count)
  colnames(count_scNPF) = colnames(count)
  count_scNPF = as.matrix(count_scNPF)
  count_scNPF
}

#Run SCRABBLE
SCRABBLE.EnImpute2=function(count,parameter = c(1,1e-6,1e-4)){
  input_count = vector('list', 2)
  input_count[[1]] = count
  count_SCRABBLE = SCRABBLE::scrabble(input_count, parameter = parameter)

  row.names(count_SCRABBLE) = row.names(count)
  colnames(count_SCRABBLE) = colnames(count)
  count_SCRABBLE = as.matrix(count_SCRABBLE)
  count_SCRABBLE
}

# Run scRMD
scRMD.EnImpute2 = function(count, tau = NULL, lambda = NULL, candidate = 0.05){

  # library and log normalization
  totalUMIPerCell = colSums(count)
  count_norm = log10(sweep(count, 2, totalUMIPerCell/1000000, '/')+1);

  count.t = t(count_norm)
  cutoff = quantile(count.t[count.t>0], candidate)
  # Impute using the function rmd
  count.scRMD  = scRMD::rmd(count.t, candidate = cutoff)$exprs
  count.scRMD = t(count.scRMD)

  count.scRMD = 10^count.scRMD-1

  row.names(count.scRMD) = row.names(count)
  colnames(count.scRMD) = colnames(count)
  count.scRMD = as.matrix(count.scRMD)
  count.scRMD
}

#Run scTSSR
scTSSR.EnImpute2 <- function(count,lambda1 = NULL, lambda2 = 1e+10, initA = NULL,
                             initB = NULL, percent = 0, ncores = 1, MAX_ITER = 4,
                             ABSTOL = 0.001, learning_rate = 1e-04, epochs = 100,
                             batch_size = 128, run_batch = TRUE, verbose = TRUE,
                             estimates.only = FALSE){
  count_scTSSR <- scTSSR::scTSSR(count, lambda1=lambda1, lambda2=lambda2,initA=initA,
                                 initB=initB, percent=percent, ncores=ncores,MAX_ITER=MAX_ITER,
                                 ABSTOL=ABSTOL,learning_rate=learning_rate, epochs= epochs,
                                 batch_size=batch_size,run_batch=run_batch,verbose=verbose,
                                 estimates.only=estimates.only)$estimate

  count_scTSSR_name <- rownames(count_scTSSR)
  count_name <- rownames(count)
  match_result <- match(count_scTSSR_name,count_name)
  count_scTSSR_name1 <- data.frame(count_scTSSR_name)
  count_name1<- data.frame(count_name)
  match_result1<- data.frame(match_result)
  count_name2 <- data.frame(x=1:nrow(count_name1))
  filter <- subset(count_name2,!count_name2$x %in% c(match_result1$match_result))
  count_scTSSR <- rbind(count_scTSSR,count[filter$x,])[count_name,]

  # row.names(count_scTSSR) <- row.names(count)
  # colnames(count_scTSSR) <- colnames(count)
  count_scTSSR <- as.matrix(count_scTSSR)
  count_scTSSR
}

#Run scTSSR2
scTSSR2.EnImpute2=function(count,k.gene = NULL, k.cell = NULL, W = NULL,
                           lambda = 256, percent = 0, ncores = 1, MAX.ITER = 4,
                           ABSTOL = 0.001, learning.rate = 1e-04, epochs = 100,
                           verbose = TRUE, estimates.only = TRUE){
  count_scTSSR2=scTSSR2::scTSSR2(count, k.gene = k.gene, k.cell = k.cell, W = W,
                                 lambda = lambda, percent = percent, ncores = ncores, MAX.ITER = MAX.ITER,
                                 ABSTOL = ABSTOL, learning.rate = learning.rate, epochs = epochs,
                                 verbose = verbose, estimates.only = estimates.only)

  row.names(count_scTSSR2) = row.names(count)
  colnames(count_scTSSR2) = colnames(count)
  count_scTSSR2 = as.matrix(count_scTSSR2)
  count_scTSSR2
}

#Run SDImpute
SDImpute.EnImpute2=function(count,do.nor=FALSE,do.log=TRUE,criterion = "asw",krange=c(5:15),k=5,M=15,T=0.5){
  count_SDImpute<-SDImpute::SDImpute(count,
                                     do.nor=do.nor,
                                     do.log=do.log,
                                     auto.k=TRUE,
                                     criterion = criterion,
                                     krange=krange,
                                     k=k,
                                     M=M,
                                     T=T)
  row.names(count_SDImpute) = row.names(count)
  colnames(count_SDImpute) = colnames(count)

  count_SDImpute = as.matrix(count_SDImpute)
  count_SDImpute
}

#Run VIPER
VIPER.EnImpute2=function(count, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1, report = FALSE, outdir = NULL, prefix = NULL){
  library(VIPER)
  count_VIPER = VIPER::VIPER(count, num = num, percentage.cutoff = percentage.cutoff, minbool = minbool, alpha = alpha,report = report, outdir = outdir, prefix = prefix)$imputed

  row.names(count_VIPER) = row.names(count)
  colnames(count_VIPER) = colnames(count)
  count_VIPER = as.matrix(count_VIPER)
  count_VIPER
}

#Run zinbwave
zinbwave.EnImpute2=function(count,nb.repeat.initialize = 2,maxiter.optimize = 25,
                      stop.epsilon.optimize = 1e-04){
  count=round(count)
  count_zinbFit <- zinbwave::zinbFit(count,nb.repeat.initialize = nb.repeat.initialize,maxiter.optimize = maxiter.optimize,
                                     stop.epsilon.optimize = stop.epsilon.optimize)
  count_zinbwave = t(zinbwave::imputeZeros(count_zinbFit, t(count)))

  row.names(count_zinbwave) = row.names(count)
  colnames(count_zinbwave) = colnames(count)
  count_zinbwave = as.matrix(count_zinbwave)
  count_zinbwave
}

#' Run EnImpute2 on a raw read count matrix
#'
#' This function is implemented to perform EnImpute2 on a raw read count matrix. EnImpute2 is an ensemble
#' learning-based method for imputing dropout values in scRNA-seq data which is a refinement of EnImpute which
#' integrates eighteen state-of-the-art methods: ALRA,DCA,DrImpute,knn_smooth, MAGIC,SAVER,
#' scImpute, scNPF,SCRABBLE,scRMD,scTSSR,scTSSR2,SDImpute,VIPER and zinbwave. EnImpute2 first runs the
#' fifteen individual imputation methods, and then use the trimmed mean of the imputed values generated
#' by different individual methods as a gold standard result.Then EnImpute2 converts all imputation matrixs into different vectors by columns
#' and calculate the Pearson correlation coefficients between eighteen methods’ imputation results and the gold standard result.
#' Finally EnImpute2 utilizes bcp package to sort out the most suitable methods for the data and takes the weighted average
#' (regarding Pearson correlation coefficients as weights)of them as the final result.This function depends on the follwing R package:
#' DrImpute, Rmagic, SAVER, scImpute, WGCNA,scNPF,SCRABBLE,scRMD,scTSSR,scTSSR2,SDImpute,VIPER,zinbwave,rsvd,bcp. These packages will be automatically installed along
#' with EnImpute2. EnImpute2 also depends on the following seven Python packages: DCA,MAGIC,scTSSR,scTSSR2. Before
#' installing the R package EnImpute2, please install the seven Python packages following the corresponding
#' readme files, and check whether they can be run from the command line.
#'
#' @param count raw read count matrix. The rows correspond to genes and the columns correspond to cells.
#' @param scale.factor scale factor used to re-scale the imputed results generated by
#' different individual methods. Default is 10000.
#' @param trim specifies the fraction (between 0 and 0.5)  of observations to be trimmed
#' from each end before the mean is computed. Default is 0.3.
#' @param ALRA a boolean variable that defines whether to impute the raw data using the ALRA method.
#' Default is TRUE.
#' @param DCA a boolean variable that defines whether to impute the raw data using the DCA method.
#' Default is TRUE.
#' @param DrImpute a boolean variable that defines whether to impute the raw data using the DrImpute method.
#' Default is "TRUE".
#' @param knn_smooth a boolean variable that defines whether to impute the raw data using the knn_smooth method.
#' Default is "TRUE".
#' @param MAGIC a boolean variable that defines whether to impute the raw data using the MAGIC method.
#' Default is TRUE.
#' @param SAVER a boolean variable that defines whether to impute the raw data using the SAVER method.
#' Default is TRUE.
#' @param scImpute a boolean variable that defines whether to impute the raw data using the scImpute method.
#' Default is TRUE.
#' @param scNPF a boolean variable that defines whether to impute the raw data using the scNPF method.
#' Default is TRUE.
#' @param SCRABBLE a boolean variable that defines whether to impute the raw data using the SCRABBLE method.
#' Default is TRUE.
#' @param scRMD a boolean variable that defines whether to impute the raw data using the scRMD method.
#' Default is TRUE.
#' @param scTSSR a boolean variable that defines whether to impute the raw data using the scTSSR method.
#' Default is TRUE.
#' @param scTSSR2 a boolean variable that defines whether to impute the raw data using the scTSSR2 method.
#' Default is TRUE.
#' @param SDImpute a boolean variable that defines whether to impute the raw data using the SDImpute method.
#' Default is TRUE.
#' @param VIPER a boolean variable that defines whether to impute the raw data using the VIPER method.
#' Default is TRUE.
#' @param zinbwave a boolean variable that defines whether to impute the raw data using the zinbwave method.
#' Default is TRUE.
#' @param ALRA.k the rank of the rank-k approximation in ALRA. Set to 0 for automated choice of k.
#' Default is 0.
#' @param ALRA.q the number of power iterations in randomized SVD used by ALRA. Default is 10.
#' @param DCA.normtype a string variable specifying the type of size factor estimation in DCA.
#' Possible values: "deseq", "zheng". Default is "zheng".
#' @param DCA.type a string variable specifying type of autoencoder in DCA. Possible values:
#' "normal", "poisson", "nb", "nb-shared", "nb-conddisp", "nb-fork", "zinb", "zinb-shared", "zinb-conddisp",
#' "zinb-fork". Default is "zinb-conddisp".
#' @param DCA.l2 a real number specifying the L2 regularization coefficient in DCA.  Default is 0.
#' @param DCA.l1 a real number specifying the L1 regularization coefficient in DCA.  Default is 0.
#' @param DCA.l2enc a real number specifying the encoder-specific L2 regularization coefficient in DCA.
#' Default is 0.
#' @param DCA.l1enc a real number specifying the encoder-specific L1 regularization coefficient in DCA.
#' Default is 0.
#' @param DCA.ridge a real number specifying the L2 regularization coefficient for dropout probabilities
#' in DCA. Default is 0.
#' @param DCA.gradclip a real number specifying the Clip grad values in DCA. Default is 5.
#' @param DCA.activation a string value specifying the activation function of hidden unit in DCA. Default is "relu".
#' @param DCA.hiddensize a string vector specifying the size of hidden layers in DCA. Default is "64,32,64".
#' @param DCA.hyper a logical value specifying whether hyperparameter search is performed in DCA.
#' @param DCA.hypern an integer specifying the number of samples drawn from hyperparameter distributions
#' during optimization in DCA. Default is 1000.
#' @param DrImpute.ks an integer vector specifying the number of cell clustering groups in DrImpute.
#' Default is 10:15.
#' @param DrImpute.dists a string vector specifying the distance metrics in DrImpute. Default is
#' c("spearman", "pearson").
#' @param DrImpute.method a string specifying the method used for imputation in DrImpute. Use "mean"
#' for mean imputation, "med" for median imputation.
#' @param DrImpute.cls a matrix specifying the clustering information manually provided by users in DrImpute.
#' The rows represent different clusterings, and the columns represent cells. Default is NULL,
#' which means the user do not provide the clustering information.
#' @param MAGIC.k an integer specifying the number of nearest neighbors on which to build kernel in MAGIC.
#' Default is 10.
#' @param MAGIC.alpha an integer specifying the decay rate of kernel tails in MAGIC. Default is 15.
#' @param MAGIC.t an integer specifying the diffusion time for the Markov Affinity Matrix in MAGIC.
#' Default is "auto". For detail about the approach to set paramter t automatically,
#' please refer to the reference.
#' @param MAGIC.npca  an integer specifying the number of PCA components in MAGIC.
#' Default is 20.
#' @param MAGIC.t.max an integer specifying the maximum value of t to test for automatic t selection in MAGIC.
#' Default is 20.
#' @param MAGIC.knn.dist.method  a string value specifying the metric for building kNN graph in MAGIC.
#' Recommended values: "euclidean", "cosine". Default is "euclidean".
#' @param MAGIC.n.jobs an integer specifying the number of jobs used for computation in MAGIC. If -1 all CPUs are used.
#' If 1 is given, no parallel computing code is used at all. For n.jobs below -1, (n.cpus + 1 + n.jobs)
#' are used. Thus for n.jobs = -2, all CPUs but one are used.
#' @param SAVER.do.fast a boolean variable specifying whether the prediction step is
#' approximated in SAVER. Default is TRUE.
#' @param SAVER.ncores number of cores to use in SAVER. Default is 1.
#' @param SAVER.size.factor a vector of cell size specifying the normalization factors in SAVER.
#' If the data is already normalized or normalization is not desired, set size.factor = 1.
#' Default uses mean library size normalization.
#' @param SAVER.npred number of genes for regression prediction in SAVER. Selects the top npred genes in
#' terms of mean expression for regression prediction. Default is all genes.
#' @param SAVER.null.model a boolean variable specifying whether to use mean gene expression as prediction
#' in SAVER. Default is FALSE
#' @param SAVER.mu matrix of prior means in SAVER.
#' @param scImpute.drop_thre  a number (between 0 and 1) specifying the threshold on dropout probability in scImpute.
#' Default is 0.5.
#' @param scImpute.Kcluster an integer specifying the number of cell subpopulations in scImpute. Default is 10.
#' @param scImpute.labeled  a boolean variable indicating whether cell type information is given in scImpute. Default is FALSE.
#' @param scImpute.labels  a character vector specifying the cell type in scImpute. Only needed when \code{labeled = TRUE}.
#' Default is NULL
#' @param scImpute.genelen an integer vector giving the length of each gene in scImpute.  Default is NULL.
#' @param scImpute.ncores an integer specifying the number of cores used for parallel computation in scImpute. Default is 1.
#' @param scNPF.network A adjacency matrix contation gene-gene interaction network. User can use priori mode or context mode. For priori mode, users can use publicly available molecular networks. In this package, we provided three human gene-gene interaction networks, including String, HumanNet and an integrated network. For context mode (default), a context-specific gene-gene network is constructed from the scRNA-seq data by WGCNA package.
#' @param scNPF.gamma A number between 0 and 1 (default: 0.5). gamma is the trade-off between prior information and network diffusion, governing the distance that a signal is allowed to diffuse through the network during smoothing. The specific value of gamma has little effect on the results of network propagation over a sizable range.
#' @param scNPF.qt.gene A numeric value between 0 and 1 (default: 0.4) indicating the top percent of expressed genes to be reserved for buliding a context-specific gene-gene network. Used only if network = "context".
#' @param scNPF.qt.cell A numeric value between 0 and 1 (default: 0.5) indicating the top percent of expressed cells to be reserved for buliding a context-specific gene-gene network. Used only if network = "context".
#' @param scNPF.nThreads The number of cores to use. Default is 1.
#' @param SCRABBLE.parameter the vector of parameters. The first parameter is the value of alpha in the mathematical model , the second one is the value of beta in the mathematical model.
#' @param scRMD.tau a non-negative real number specifying the tuning parameter to penalize the sparse term. Default is NULL.
#' @param scRMD.lambda a non-negative real number specifying the tuning parameter to penalize the row rank term. Default is NULL.
#' @param scRMD.candidate a real number (0 to 1) specifying the cutoff for candidate drop out. Default is 0.05.
#' @param scTSSR.lambda1 Tuning parameter to facilitate feature selection and regularization.
#' @param scTSSR.lambda2 Tuning parameter to penalize the diagonal elements of the parameter to eliminate the trivial solution of representing an expression level as a linear combination of itself.
#' @param scTSSR.initA The initionlization of A. The elements of A represent the similarities between genes.
#' @param scTSSR.initB The initionlization of B. The elements of B represent the similarities between cells.
#' @param scTSSR.percent The expression count matrix is preprocessed by filtering out the genes expressed in at most percent*100\% of the cells.
#' @param scTSSR.MAX_ITER Maximum iteration of the external circulation of scTSSR.
#' @param scTSSR.ABSTOL Absolute tolerance of the external circulation.
#' @param scTSSR.learning_rate A hyper-parameter that controls the speed of adjusting the weights of the network with respect to the loss gradient.
#' @param scTSSR.epochs The number of the entire training set going through the entire network.
#' @param scTSSR.batch_size The number of examples that are fed to the algorithm at a time.
#' @param scTSSR.run_batch Whether to use batch or to set the number of all the samples as the value of the batch size. Default is TRUE.
#' @param scTSSR.verbose Whether to output the value of metrics at the end of each epoch. Default is TRUE.
#' @param scTSSR2.k.gene A hyper-parameter that controls the sparsity level of the estimated coefficient matrices, A1 and A2. Default is k_gene = min(100, m/30).
#' @param scTSSR2.k.cell A hyper-parameter that controls the sparsity level of the estimated coefficient matrices, B1 and B2. Default is k_cell = min(100, n/30).
#' @param scTSSR2.W A weight matrix with element W_gc denotes the non-dropout probability of the expression level of gene g in cell c. Default is W_gc=X_gc/max(X_gc).
#' @param scTSSR2.lambda Ridge penalty parameter. Default is 256.
#' @param scTSSR2.percent The expression count matrix is preprocessed by filtering out the genes expressed in at most percent*100\% of the cells. Default is 0.05.
#' @param scTSSR2.ncores Number of cores to use. Default is 1.
#' @param scTSSR2.MAX.ITER Maximum iteration of the external circulation of scTSSR2. Default is 4.
#' @param scTSSR2.ABSTOL Absolute tolerance of the external circulation. Default is 1e-3.
#' @param scTSSR2.learning.rate A hyper-parameter that controls the speed of adjusting the weights of the network with respect to the loss gradient. Default is 0.0001.
#' @param scTSSR2.epochs The number of the entire training set going through the entire network. Default is 100.
#' @param scTSSR2.verbose Whether to output the value of metrics at the end of each epoch. Default is TRUE.
#' @param SDImpute.do.nor Logical. If TRUE, the data is Normalized.
#' @param SDImpute.do.log Logical. If TRUE, the input will take log of data.
#' @param SDImpute.criterion One of "asw" or "ch". Determines whether average silhouette width or Calinski-Harabasz is applied.
#' @param SDImpute.krange Integer vector. Numbers of clusters which are to be compared by the average silhouette width criterion. Note: average silhouette width and Calinski-Harabasz can't estimate number of clusters nc=1.
#' @param SDImpute.k Integer. The number of cell clusters. This parameter can be determined based on prior knowledge or clustering result of raw data.
#' @param SDImpute.M Integer. The number of nearest neighbors.When the number of nearest neighbors for each cell is small, the parameter M should not be too large to guarantee that it makes sense. In general, this parameter is set to an integer between 10 and 30.
#' @param SDImpute.T Numeric between 0 and 1. The dropout probability candidate threshold which controls the degree of imputation to the gene expression matrix. The recommended value of parameter T is 0.5.
#' @param VIPER.num The number of random sampled genes used to fit the penalized regression model to identify the set of candidate cells. The default value is 5000. If gene number p in the dataset is less than specified num, numwill be set as 0.8*p.
#' @param VIPER.percentage.cutoff To reduce the influence of missing values in the weight estimation, the nonnegative regression model is fitted using genes with a zero rate less than a certain threshold. The default value is 0.1 (10 percent).
#' @param VIPER.minbool The criteria used to select the penalty levellambda in the penalized regression model. VIPER calls cv.glmnet() in glmnet to perform fitting cross validation. Two penalty levels are available for selection: lambda.min, the value of lambda that gives minimum mean cross-validated error, and lambda.1se, the value of lambda that gives the most regularized model such that error is within one standard error of the minimum. The default is lambda.1se, i.e., minbool = FALSE.
#' @param VIPER.alpha The elastic net mixing parameter. The default value is 1, which is equivalent to a lasso model.
#' @param zinbwave.nb.repeat.initialize Number of iterations for the initialization of beta_mu and gamma_mu.
#' @param zinbwave.maxiter.optimize maximum number of iterations for the optimization step (default 25).
#' @param zinbwave.stop.epsilon.optimize stopping criterion in the optimization step, when the relative gain in likelihood is below epsilon (default 0.0001).
#'
#' @return a list with the following components
#' \item{\code{count.EnImpute2.log}}{Imputed count matrix generated by EnImpute2 (log scale).}
#' \item{\code{count.EnImpute2.exp}}{Imputed count matrix generated by EnImpute2 (exp scale).}
#' \item{\code{count.imputed.individual.exp}}{Imputed count matrices generated by different individual imputation methods (exp scale).}
#' \item{\code{Methods.used}}{The individual methods used by EnImpute.}
#'
#' @export
#' @import DrImpute Rmagic SAVER scImpute WGCNA scNPF SCRABBLE scRMD scTSSR scTSSR2 SDImpute VIPER zinbwave rsvd bcp
#' @author Yu-Heng Luo  <luoyuheng@mail.ccnu.edu.cn>
#' @examples
#'
#' data("baron")
#' result <- EnImpute2(baron$count.samp)
EnImpute2 = function(count, scale.factor = 10000,trim=0.3,threshold=0.7,
                     ALRA = TRUE, DCA=TRUE,DrImpute = TRUE,knn_smooth=TRUE,MAGIC = TRUE,
                     SAVER=TRUE,scImpute = TRUE,scNPF=TRUE,SCRABBLE=TRUE,scRMD = TRUE,scTSSR=TRUE,scTSSR2=TRUE,SDImpute=TRUE,
                     VIPER=TRUE,zinbwave=TRUE,
                     ALRA.k = 0, ALRA.q = 10,
                     DCA.normtype = "zheng", DCA.type = "zinb-conddisp",
                     DCA.l2 = 0, DCA.l1 =0, DCA.l2enc = 0, DCA.l1enc = 0, DCA.ridge = 0,
                     DCA.gradclip = 5, DCA.activation = "relu", DCA.hiddensize = "64,32,64",
                     DCA.hyper = FALSE, DCA.hypern = 1000,
                     DrImpute.ks = 10:15, DrImpute.dists = c("spearman", "pearson"),DrImpute.method = "mean", DrImpute.cls = NULL,
                     knn_smooth.k=10,
                     MAGIC.k = 10, MAGIC.alpha = 15, MAGIC.t = "auto", MAGIC.npca = 20,
                     MAGIC.t.max = 20, MAGIC.knn.dist.method = "euclidean", MAGIC.n.jobs = 1,
                     SAVER.do.fast = TRUE, SAVER.ncores = 1, SAVER.size.factor = NULL,
                     SAVER.npred = NULL, SAVER.null.model = FALSE, SAVER.mu = NULL,
                     scImpute.drop_thre = 0.5, scImpute.Kcluster = 5, scImpute.labeled = FALSE,
                     scImpute.labels = NULL, scImpute.genelen = NULL, scImpute.ncores = 1,
                     scNPF.network="context",scNPF.gamma = 0.5, scNPF.qt.gene = 0.4, scNPF.qt.cell = 0.5,
                     scNPF.nThreads = 1,
                     SCRABBLE.parameter=c(1,1e-6,1e-4),
                     scRMD.tau = NULL, scRMD.lambda = NULL, scRMD.candidate = 0.05,
                     scTSSR.lambda1 = NULL, scTSSR.lambda2 = 1e+10, scTSSR.initA = NULL,
                     scTSSR.initB = NULL, scTSSR.percent = 0.05, scTSSR.MAX_ITER = 4,
                     scTSSR.ABSTOL = 0.001, scTSSR.learning_rate = 1e-04, scTSSR.epochs = 100,
                     scTSSR.batch_size = 128, scTSSR.run_batch = TRUE, scTSSR.verbose = TRUE,scTSSR.estimates.only = FALSE,
                     scTSSR2.k.gene = NULL, scTSSR2.k.cell = NULL, scTSSR2.W = NULL,
                     scTSSR2.lambda = 256, scTSSR2.percent = 0, scTSSR2.ncores = 1, scTSSR2.MAX.ITER = 4,
                     scTSSR2.ABSTOL = 0.001, scTSSR2.learning.rate = 1e-04, scTSSR2.epochs = 100,
                     scTSSR2.verbose = TRUE, scTSSR2.estimates.only = TRUE,
                     SDImpute.do.nor=FALSE,SDImpute.do.log=TRUE,SDImpute.criterion = "asw",SDImpute.krange=c(5:15),SDImpute.k=5,SDImpute.M=15,SDImpute.T=0.5,
                     VIPER.num=5000,VIPER.percentage.cutoff=0.1, VIPER.minbool = FALSE,VIPER.alpha = 0.5,
                     zinbwave.nb.repeat.initialize = 2,zinbwave.maxiter.optimize = 25,zinbwave.stop.epsilon.optimize = 1e-04){
  Methods = c("ALRA","DCA","DrImpute","knn_smooth","MAGIC","SAVER","scImpute","scNPF","SCRABBLE","scRMD","scTSSR","scTSSR2","SDImpute","VIPER","zinbwave")
  Methods.idx = c(ALRA,DCA,DrImpute,knn_smooth,MAGIC,SAVER,scImpute,scNPF,SCRABBLE,scRMD,scTSSR,scTSSR2,SDImpute,VIPER,zinbwave)

  if(sum(Methods.idx)==0)
    stop("You need choose at least one individual imputation method.")
  Methods.used = Methods[Methods.idx]

  K = length(Methods.used)
  p = dim(count)[1]
  n = dim(count)[2]

  count.imputed.individual = array(0, dim=c(p,n,K))
  dimnames(count.imputed.individual)[[1]] = rownames(count)
  dimnames(count.imputed.individual)[[2]] = colnames(count)
  dimnames(count.imputed.individual)[[3]]= Methods.used

  k = 1
  #  ALRA
  if (ALRA == TRUE){
    count.imputed.individual[,,k]  = ALRA.EnImpute2(count, k = ALRA.k, q = ALRA.q)
    k = k +1
    print("ALRA is done!")
  }
  # DCA
  if (DCA == TRUE){
    count.imputed.individual[,,k] = DCA.EnImpute2(count, normtype = DCA.normtype, type = DCA.type,
                                                  l2 = DCA.l2, l1 = DCA.l1, l2enc = DCA.l2enc, l1enc = DCA.l1enc, ridge = DCA.ridge,
                                                  gradclip = DCA.gradclip, activation = DCA.activation, hiddensize = DCA.hiddensize,
                                                  hyper = DCA.hyper, hypern = DCA.hypern)
    k = k +1
    print("DCA is done!")
  }
  #  DrImpute
  if (DrImpute == TRUE){
    count.imputed.individual[,,k] = DrImpute.EnImpute2(count, ks = DrImpute.ks, dists = DrImpute.dists, method = DrImpute.method,
                                                       cls = DrImpute.cls)
    k = k +1
    print("DrImpute is done!")
  }
  # knn_smooth
  if (knn_smooth == TRUE){
    count.imputed.individual[,,k]  = knnsmooth.EnImpute2(count,k=knn_smooth.k)
    k = k +1
    print("knn_smooth is done!")
  }
  #  MAGIC
  if (MAGIC == TRUE){
    count.imputed.individual[,,k] = MAGIC.EnImpute2(count, k = MAGIC.k, alpha = MAGIC.alpha, t = MAGIC.t,
                                                    npca = MAGIC.npca, t.max = MAGIC.t.max,
                                                    knn.dist.method = MAGIC.knn.dist.method, n.jobs = MAGIC.n.jobs)
    k = k +1
    print("MAGIC is done!")
  }
  #  SAVER
  if (SAVER == TRUE){
    count.imputed.individual[,,k] = SAVER.EnImpute2(count, do.fast = SAVER.do.fast, ncores = SAVER.ncores,
                                                    size.factor = SAVER.size.factor, npred = SAVER.npred,
                                                    null.model = SAVER.null.model, mu = SAVER.mu)
    k = k +1
    print("SAVER is done!")
  }
  #  scImpute
  if (scImpute == TRUE){
    count.imputed.individual[,,k] = scImpute.EnImpute2(count, drop_thre = scImpute.drop_thre, Kcluster = scImpute.Kcluster,
                                                       labeled = scImpute.labeled, labels = scImpute.labels,
                                                       genelen = scImpute.genelen, ncores = scImpute.ncores)
    k = k +1
    print("scImpute is done!")
  }
  # scNPF
  if (scNPF == TRUE){
    count.imputed.individual[,,k] = scNPF.EnImpute2(count,network=scNPF.network,gamma =scNPF.gamma, qt.gene =scNPF.qt.gene,qt.cell = scNPF.qt.cell,
                                                    nThreads =scNPF.nThreads)
    k = k +1
    print("scNPF is done!")
  }
  # SCRABBLE
  if (SCRABBLE == TRUE){
    count.imputed.individual[,,k] = SCRABBLE.EnImpute2(count,parameter = SCRABBLE.parameter)
    k = k +1
    print("SCRABBLE is done!")
  }
  #  scRMD
  if (scRMD == TRUE){
    count.imputed.individual[,,k] = scRMD.EnImpute2(count, tau = scRMD.tau, lambda = scRMD.lambda,
                                                    candidate = scRMD.candidate)
    k = k +1
    print("scRMD is done!")
  }
  #  scTSSR
  if (scTSSR == TRUE){
    count.imputed.individual[,,k] = scTSSR.EnImpute2(count,lambda1 = scTSSR.lambda1, lambda2 = scTSSR.lambda2, initA = scTSSR.initA,
                                                     initB = scTSSR.initB, percent = scTSSR.percent, ncores = 1, MAX_ITER = scTSSR.MAX_ITER,
                                                     ABSTOL = scTSSR.ABSTOL, learning_rate = scTSSR.learning_rate, epochs = scTSSR.epochs,
                                                     batch_size = scTSSR.batch_size, run_batch = scTSSR.run_batch, verbose = scTSSR.verbose,
                                                     estimates.only = scTSSR.estimates.only)
    k = k +1
    print("scTSSR is done!")
  }
  # scTSSR2
  if (scTSSR2 == TRUE){
    count.imputed.individual[,,k] = scTSSR2.EnImpute2(count,k.gene = scTSSR2.k.gene, k.cell = scTSSR2.k.cell, W = scTSSR2.W,
                                                      lambda = scTSSR2.lambda, percent = scTSSR2.percent, ncores = scTSSR2.ncores, MAX.ITER = scTSSR2.MAX.ITER,
                                                      ABSTOL = scTSSR2.ABSTOL, learning.rate = scTSSR2.learning.rate, epochs = scTSSR2.epochs,
                                                      verbose = scTSSR2.verbose, estimates.only = scTSSR2.estimates.only)
    k = k +1
    print("scTSSR2 is done!")
  }
  #  SDImpute
  if (SDImpute == TRUE){
    count.imputed.individual[,,k] = SDImpute.EnImpute2(count,do.nor=SDImpute.do.nor,do.log=SDImpute.do.log,criterion = SDImpute.criterion,krange=SDImpute.krange,k=SDImpute.k,M=SDImpute.M,T=SDImpute.T)
    k = k +1
    print("SDImpute is done!")
  }
  #  VIPER
  if (VIPER == TRUE){
    count.imputed.individual[,,k] = VIPER.EnImpute2(count,num=VIPER.num,percentage.cutoff=VIPER.percentage.cutoff, minbool = VIPER.minbool,alpha = VIPER.alpha)
    k = k +1
    print("VIPER is done!")
  }
  #  zinbwave
  if (zinbwave == TRUE){
    count.imputed.individual[,,k] = zinbwave.EnImpute2(count,nb.repeat.initialize = zinbwave.nb.repeat.initialize,maxiter.optimize = zinbwave.maxiter.optimize,
                                                       stop.epsilon.optimize = zinbwave.stop.epsilon.optimize)
    k = k +1
    print("zinbwave is done!")
  }
  k=k-1

  count.imputed.individual[count.imputed.individual<=0] = 0
  count.imputed.individual[is.na(count.imputed.individual)] = 0

  count.imputed.individual.rescaled = count.imputed.individual # Rescale the imputed count matrices
  for (i in 1:k){
    totalUMIPerCell = colSums(count.imputed.individual[,,i])
    count.imputed.individual.rescaled[,,i] = sweep(count.imputed.individual[,,i] , 2, totalUMIPerCell/scale.factor, '/');
  }
  count.imputed.individual.rescaled = log(count.imputed.individual.rescaled+1)

  count.EnImpute2.log = apply(count.imputed.individual.rescaled, 1:2, mean,trim=trim)
  rownames(count.EnImpute2.log) = rownames(count)
  colnames(count.EnImpute2.log) = colnames(count)
  count.EnImpute2.exp = exp(count.EnImpute2.log) - 1#obtain the standard result

  first_vec=as.vector(count.EnImpute2.exp)#correlation calculation
  all_cor=matrix(0,1,k)
  count_vec <-matrix(0,k,length(count.EnImpute2.exp))
  for (i in 1:k){
    count_vec[i,] <-as.vector(count.imputed.individual[,,i])
  }
  for(i in 1:k){
    all_cor[,i]=cor(count_vec[i,],first_vec)
  }
  all_cor=all_cor/sum(all_cor)#归一化

  new=data.frame(index=c(1:k),cor=all_cor[1,])#get the best k
  new1=new[order(new$cor,decreasing = TRUE),]
  new1$index=c(1:k)

  bcp_x <- bcp::bcp(new1$cor, return.mcmc = TRUE)#bcp change point detection
  bcp_sum <- as.data.frame(summary(bcp_x))
  bcp_sum$id <- 1:length(new1$cor)
  sel <- bcp_sum[which(bcp_x$posterior.prob > threshold), ]
  if(is.na(sel$id[1])){
    change_point=k
  }else{
    change_point=sel$id[1]
  }
  method_ind=t(apply(all_cor,1,order,decreasing=TRUE)[1:change_point,])

  count.EnImpute2.log = apply(count.imputed.individual.rescaled[,,c(method_ind)], 1:2, weighted.mean,w=all_cor[c(method_ind)])#obtain the final result

  rownames(count.EnImpute2.log) = rownames(count)
  colnames(count.EnImpute2.log) = colnames(count)
  count.EnImpute2.exp = exp(count.EnImpute2.log) - 1
  result = list(count.EnImpute2.log = count.EnImpute2.log, count.EnImpute2.exp = count.EnImpute2.exp,
                count.imputed.individual = count.imputed.individual, count.imputed.individual.rescaled = count.imputed.individual.rescaled,
                Methods.used = Methods.used)
  result
}




