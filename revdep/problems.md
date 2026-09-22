# iCAMP (1.8.6)

* GitHub: <https://github.com/DaliangNing/iCAMP1>
* Email: <mailto:ningdaliang@ou.edu>
* GitHub mirror: <https://github.com/cran/iCAMP>

Run `revdepcheck::revdep_details(, "iCAMP")` for more info

## In both

*   checking examples ... ERROR
     ```
     ...
     > 
     > data("example.data")
     > comm=example.data$comm
     > pd=example.data$pd
     > 
     > # in this example, 10 samples from one metacommunity,
     > # the other 10 samples from another metacommunity.
     > meta.group=data.frame(meta.com=c(rep("meta1",10),rep("meta2",10)))
     > rownames(meta.group)=rownames(comm)
     > 
     > nworker=2 # parallel computing thread number.
     > rand.time=20 # usually use 1000 for real data.
     > sigmpd=NRI.cm(comm=comm, meta.group=meta.group,
     +               dis=pd, nworker=nworker,
     +               weighted=TRUE, rand=rand.time,
     +               sig.index="all")
     Warning: 'memory.limit()' is no longer supported
     All match very well.
     All match very well.
     Now calculating observed MPD. Begin at Tue Sep 22 15:24:16 2026. Please wait...
     Now randomizing by parallel computing. Begin at Tue Sep 22 15:24:17 2026. Please wait...
     Error in checkForRemoteErrors(val) : 
       2 nodes produced errors; first error: there is no package called 'iCAMP'
     Calls: NRI.cm ... clusterApply -> staticClusterApply -> checkForRemoteErrors
     Execution halted
     ```

