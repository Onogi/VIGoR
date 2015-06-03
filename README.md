# VIGoR
VIGoR (Variational Bayesian Inference for Genome-Wide Regression) conducts linear regression using variational Bayesian inference, particularly optimized for genome-wide association mapping and whole-genome prediction which use a number of DNA markers as the explanatory variables. VIGoR provides seven regression models which select the important variables (i.e., the variables related to response variables) among the given explanatory variables in different ways (i.e., model structures). This directry contains the followings.

VIGoR_document.pdf: Document file for VIGoR.

VIGoR_Linux: Command line program package for Linux. This was compiled under the Linux kernel release 3.13.0-24-generic with a X86-64   machine.

VIGoR_Mac: Command line program package for Mac. This was compiled under OS X ver. 10.6.8.

Source: C source files of the command line programs.

R: R source package of VIGoR. Installation of R source packages will be easy on Linux and Mac. However, it is often problematic on Windows. I recommend the following procedure for the installation on Windows:
  1. Install Rtools from http://cran.r-project.org/bin/windows/Rtools/
  2. Install Rstudio from http://www.rstudio.com/
  3. Start up Rstudio and install the R source package (.tar.gz file) as a package archive file. The package compiled by Rstudio can be used from the standard R console.

Copyright (C) 2015 Akio Onogi and Hiroyoshi Iwata.
Released under the MIT license.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

