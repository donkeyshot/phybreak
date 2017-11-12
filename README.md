# phybreak
Outbreak reconstruction with sequence data

The package implements the method described in Klinkenberg et al (2016), doi: http://dx.doi.org/10.1101/069195

Workflow:

* enter data and priors by constructing an object of S3-class 'phybreak', with function 'phybreak'

* do mcmc-updates with functions 'burnin_phybreak' and 'sample_phybreak'; remove samples with 'thin.phybreak'

* access the 'phybreak'-object by get_phybreak-functions such as 'get_transtree', 'get_data', 'get_parameters'

* summarize the mcmc-chain with the functions 'ESS', 'transtree', 'infectorsets', 'phylotree', 'get_mcmc', 'get_phylo'

* plotting with 'plot', 'plotTrans', and 'plotPhylo'


* it is possible to simulate data with 'sim_phybreak'