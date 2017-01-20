# phybreak 0.1.1.9004

## BUG FIXES

* Bug in axis plotting with plotTrans if time did not start at 0
* Bug in sim.phybreak giving error with option output.class = "obkData"


# phybreak 0.1.1.9003

## Minor Changes

* Faster simulation and phybreak construction with long sequences.


# phybreak 0.1.1.9002

## Major Changes

* Added mcmc-update steps for only the within-host phylogenetic minitree topology (not used by default)

## Minor Changes

* Better user control of plotTrans graphical parameters


# phybreak 0.1.1.9001

## Major Changes

* Added transmission tree plotting function plotTrans
* Data entry through new S3-class phybreakdata
* Simulation with sim.phybreak results in phybreakdata-object, that can be plotted

## Minor Changes

* Automatic calculation of initial value for mu
* Reference date in phybreak object (d-slot) referring to t = 0
* Minor changes in get.phybreak functions


# phybreak 0.1.1

## BUG FIXES

* Solves issue with irreconcilable MCC trees in phylo-class vs table-form (#1) (reported by Mark Schultz)  
* Solves errors with OS X 10.9 Mavericks

