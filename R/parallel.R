#################GENERIC PARALLEL INTERFACE #################

DBA_PARALLEL_MULTICORE = 1
DBA_PARALLEL_RLSF      = 2

### INITIALIZE ###
dba.parallel = function(DBA) {

   if(DBA$config$parallelPackage==DBA_PARALLEL_MULTICORE) {
      DBA$config =  dba.multicore.init(DBA$config)
      DBA$config$parallelInit = TRUE
      return(DBA)
   }
   
   if(DBA$config$parallelPackage==DBA_PARALLEL_RLSF) {
      DBA$config =  dba.Rlsf.init(DBA$config)
      DBA$config$parallelInit = TRUE
      return(DBA)
   }
   
   warning('UNKNOWN PARALLEL PACKAGE: ',DBA$config$parallelPackage)
   return(DBA)
}

dba.parallel.init = function(config){
   config = config$initFun(config)
   return(config)
}

### CONFIG ###
dba.parallel.params = function(config,funlist) {
   params = config$paramFun(config,funlist)
   return(params)
}

### PARALLEL LAPPLY ###
dba.parallel.lapply = function(config,params,arglist,fn,...) {
   jobs = config$lapplyFun(config,params,arglist,fn,...)
   return(jobs)
}

### ADD JOB ###
dba.parallel.addjob = function(config,params,fn,...) {
   job = config$addjobFun(config,params,fn,...)
   return(job)
}

### WAIT FOR JOBS TO COMPLETE AND GATHER RESULTS ###
dba.parallel.wait4jobs = function(config,joblist) {
   res = config$wait4jobsFun(config,joblist)
   return(res)
}

dba.Rlsf.init = function(config){
   if (length(find.package(package="DiffBindCRI",quiet=T))>0) {
      #library(DiffBindCRI)
      #if(!exists("dba.Rlsf.init","package:DiffBindCRI")) {
      #   warning('Rlsf interface not supported in this version')
      #} else {
         config = dba.CRI.Rlsf.init(config)
      #}     
  } else {
    warning('Rlsf interface not supported in this version')
  }
   return(config)
}

################# multicore INTERFACE #################

### INITIALIZE ###

dba.multicore.init = function(config) {

   noparallel=F
   if (length(find.package(package="parallel",quiet=T))>0) {
      library(parallel)
      if(!exists("mcparallel","package:parallel")) {
         noparallel=T
      }     
      if(!exists("mclapply","package:parallel")) {
         noparallel=T
      } 
   } else {
      noparallel=T
   }
   if(noparallel){
      warning("Parallel execution unavailable: executing serially.")
      config$RunParallel = FALSE
      config$parallelPackage = 0
      return(config)
   }
   
   config$multicoreInit = T    

   config$initFun      = dba.multicore.init
   config$paramFun     = dba.multicore.params
   config$addjobFun    = dba.multicore.addjob
   config$lapplyFun    = dba.multicore.lapply
   config$wait4jobsFun = dba.multicore.wait4jobs

   config$cores = parallel:::detectCores(logical=FALSE)
   
   return(config)
}

### CONFIG ###
dba.multicore.params = function(config,funlist) {
   return(NULL)
}

### PARALLEL LAPPLY ###
dba.multicore.lapply = function(config,params,arglist,fn,...){
   savecores = options("cores")
   savemccores = options("mc.cores")
   options(cores = config$cores)
   options(mc.cores = config$cores)
   res = mclapply(arglist,fn,...,mc.preschedule=TRUE) #changed to TRUE because of bug in R 2.15
   options(cores=savecores)
   options(mc.cores=savemccores)
   return(res)
}

### ADD JOB ###
dba.multicore.addjob = function(config,params,fn,...) {
  job = mcparallel(fn(...))
  return(job)
}

### WAIT FOR JOBS TO COMPLETE AND GATHER RESULTS ###
dba.multicore.wait4jobs = function(config,joblist) {
   res = mccollect(joblist)
   return(res)
}
