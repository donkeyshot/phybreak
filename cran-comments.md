### Test environments  
* Ubuntu 16.4, R 3.4.1  
* win-builder  

### R CMD check results
There were no ERRORs or WARNINGS

There was 1 NOTE:

***  
\* checking CRAN incoming feasibility ... NOTE  
Maintainer: 'Don Klinkenberg <don@xs4all.nl>'  

Version contains large components (0.2.0)  

Possibly mis-spelled words in DESCRIPTION:  
    Klinkenberg (9:5)  
    al (9:20)  
    datasets (10:33)  
    et (9:17)  

***  
This may be due to the vignette (200kb), or data (300kb)
