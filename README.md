CSSP: ChIP-SEQ Statistical Power
================================

What is CSSP?
-------------

CSSP (*ChIP-SEQ Statistical Power*) is a statistical framework for computing the power of detecting enriched regions for ChIP-SEQ experiments.

How does it Work?
-----------------

CSSP framework incorporates the following steps:

- Constructing bin-level data from aligned short reads for both a ChIP sample and a control sample;
- Preparing bin-level genome mappability and gc-content scores;
- Fitting the proposed local Poisson model for the bin-level sequencing data;
- Calculating power at desired FDR and effect sizes.

Installation
------------
CSSP is currently available in bioconductor. A stable version can be accessed through Bioconductor by:

     source( "http://cran.cnr.berkeley.edu" )
     biocLite( "CSSP" )

References
----------

Zuo, C. and Keles, S., A Statistical Framework for Power Calculations in ChIP-Seq Experiments, Bioinformatics(2013) doi: 10.1093/bioinformatics/btt200

Acknowledgement
---------------

CSSP is developed by Keles Lab at University of Wisconsin Madison. 

