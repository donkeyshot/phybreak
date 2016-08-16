# phybreak
Outbreak reconstruction with sequence data

The package implements the method described in Klinkenberg et al (2016), doi: http://dx.doi.org/10.1101/069195

Workflow:
* enter data and priors by constructing an object of S3-class 'phybreak', with function 'phybreak'
* do mcmc-updates with functions 'burnin.phybreak' and 'sample.phybreak'; remove samples with 'thin.phybreak'
* access the 'phybreak'-object by get.phybreak-functions such as 'get.tree', 'get.seqdata', 'get.parameters'
* summarize the mcmc-chain with the functions 'transtree', 'infectorsets', 'phylotree', 'get.mcmc', 'get.phylo'

* it is possible to simulate data with 'sim.phybreak'