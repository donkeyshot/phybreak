# phybreak 0.3.0

### Major changes
* Added three within-host models allowing for wide transmission bottlenecks, i.e. bottlenecks with more than 1 lineage:
    * Change in variable and sample slots of phybreak objects, with separate infector and inftime entries for the transmission tree
    * More within-host parameters in the parameter slot of phybreak objects
    * Simulation and analysis with the new models
    * See help(phybreak) for a description of the new models
* Effective sample size calculation with new function ESS
    * Applied to factor objects, it applies a method similar to the approximate ESS for phylogenetic tree (Lanfaer et al, 2016)
    * Applied to phybreak objects, it applies ESS.factor for the infectors, and coda::effectiveSize for all other parameters and variables
    * See help(ESS) for more information
* More generic functions working with phybreak objects: c, print, quantile, thin
* Dots in function names replaced by underscores
* New update protocol for proposing the within-host minitrees. In the "classic" protocol, each proposal of a new infection time and infector for a host included a complete resampling of its minitree. In the new protocol (which is now the default in "burnin_phybreak" and "sample_phybreak") only new coalescent times are proposed, but the topology is not. This is now followed by proposals in which the sampling tips are one by one removed and reattached (by simulation). 

### Minor changes
* Function infectorsets can produce matrix with support for all host-infector pairs
* New arguments in "burnin_phybreak" and "sample_phybreak":
    * classic: to choose the proportion of tree updates with the classic protocol instead of the new default protocol (see above)
    * parameter_frequency: option to sample the model parameters more frequently than each host in the transmission tree
    * status_interval: the time interval between the on-screen status updates of the mcmc-chain
* function get_bottlenecks to obtain the number of lineages leaving each host to its infector

### BUG FIXES
* Fixed bug in transtree-method "edmonds" that could occur with much uncertainty in the posterior trees


# phybreak 0.2.1

### BUG FIXES
* Fixed bug in matching samples to hosts in phybreakdata (#4 reported by Gerry Tonkin-Hill)


# phybreak 0.2.0

### Major changes
* Datasets added
* Vignette phybreak_intro added
* Use more samples per host, with changes in data input, the phybreak object structure, and plotting
* Added mcmc-update steps for only the within-host phylogenetic minitree topology (not used by default)
* Added transmission tree plotting function plotTrans
* Data entry through new S3-class phybreakdata
* Simulation with sim.phybreak results in phybreakdata-object, that can directly be plotted
* MCMC progress shown on screen (10 sec intervals)

### Minor changes
* Faster simulation and phybreak construction with long sequences.
* Automatic calculation of initial value for mu
* Reference date in phybreak object (d-slot) referring to t = 0, used for transtree and plot
* Minor changes in get.phybreak functions

### BUG FIXES
* Fixed bug in transtree with Edmonds method (index case could have been incorrect)
* Correct on-screen log-likelihood during MCMC
* Bug in sim.phybreak giving error with option output.class = "obkData"


# phybreak 0.1.1

### BUG FIXES

* Solves issue with irreconcilable MCC trees in phylo-class vs table-form (#1) (reported by Mark Schultz)  
* Solves errors with OS X 10.9 Mavericks

