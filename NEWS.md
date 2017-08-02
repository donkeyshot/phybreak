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

