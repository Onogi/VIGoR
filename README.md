# VIGoR
VIGoR (Variational Bayesian Inference for Genome-Wide Regression) conducts linear regression using variational Bayesian inference, particularly optimized for genome-wide association mapping and whole-genome prediction which use a number of DNA markers as the explanatory variables. VIGoR provides seven regression models which select the important variables (i.e., the variables related to response variables) among the given explanatory variables in different ways (i.e., model structures). This directry contains the followings.

VIGoR_document.pdf:
  Document file for VIGoR.

VIGoR_Linux:
  Command line program package for Linux. This was compiled under the Linux kernel release 3.13.0-24-generic with a X86-64   machine.

VIGoR_Mac:
  Command line program package for Mac. This was compiled under OS X ver. 10.6.8.

Source:
  C source files of the command line programs.

R:
  R source package of VIGoR. Whereas installation of R source packages is easy on Linux and Mac, it is often problematic on Windows. I recommend the following procedure:
  1. Install Rtools from http://cran.r-project.org/bin/windows/Rtools/
  2. Install Rstudio from http://www.rstudio.com/
  3. Start up Rstudio and install the R source package (.tar.gz file) as a package archive file.
  The package compiled by Rstudio can be used from the standard R console.
