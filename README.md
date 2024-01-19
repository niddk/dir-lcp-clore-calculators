# Clore Lab Protein Calculators

This repository contains source code for protein calculators developed by the G. Marius Clore lab
([Section of Molecular and Structural Biophysics][section]) of the [Laboratory of Chemical Physics][lab] at the [National Institute
of Diabetes and Digestive and Kidney Diseases (NIDDK)][niddk].

## Overview

These are the files for setup of various protein calculators for use as referenced in respective publications (**please cite when using these programs**):

- Paramagnetic Relaxation Enhancement (PRE) Calculator (simple, model free, and extended model free versions):
  - Anthis NJ, Clore GM. Visualizing transient dark states by NMR spectroscopy. _Q Rev Biophys._ 2015;48(1):35-116. doi:10.1017/S0033583514000122
  - [PubMed](https://www.pubmed.ncbi.nlm.nih.gov/25710841)
- A205 Protein/Peptide Concentration Calculator:
  - Anthis NJ, Clore GM. Sequence-specific determination of protein and peptide concentrations by absorbance at 205 nm. _Protein Sci._ 2013;22(6):851-858. doi:10.1002/pro.2253
  - [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/23526461)

## Inventory

This repository contains web content (HTML and PDF) as well as CGI applications, implemented as Python scripts.

- `/cgi-bin` - CGI applications to power the calculators
  - `AbsCalc.cgi` - A205 Protein/Peptide Concentration Calculator logic (Python)
  - `pre-calc.cgi` - PRE Calculator (Simple Version) logic (Python)
  - `pre-calc_MF.cgi` - PRE Calculator (Model Free Version) logic (Python)
  - `pre-calc_extMF.cgi` - PRE Calculator (Extended Model Free Version) logic (Python)
- `/html` - HTML pages providing a user interface for the calculators
  - `A205.html` - A205 Protein/Peptide Concentration Calculator form page
  - `precalc.html` - PRE Calculator (Simple Version) form page
  - `precalc-mf.html` - PRE Calculator (Model Free Version) form page
  - `precalc-extmf.html` - PRE Calculator (Extended Model Free Version) form page
- `/html/pdf` - Related publications, linked from the HTML pages
  - `473.pdf` - PRE Calculator (Simple, Model Free, and Extended Model Free Versions) publication
  - `453.pdf` - A205 Protein/Peptide Concentration Calculator publication

## Setup

1. Install and configure [Apache server][apache]
2. Enable the Apache CGI module
3. Install Python
4. Add the Clore CGI scripts to the `cgi-bin` directory
5. (If on Unix/Linux) Mark the `.cgi` scripts as executable (_e.g._, `chmod +x *.cgi`)
6. Add the contents of the `html` directory to the the web root (_e.g._, `/var/www/html`)

## Contact Information and Resources

- [Dr. Clore's Email][email]
- [Dr. Clore's Biography][bio]
- [Dr. Clore's Select Publications][pubs]

## License

This code is copyright free and no license is required for its use.

[email]: mailto:mariusc@intra.niddk.nih.gov
[bio]: https://www.niddk.nih.gov/about-niddk/staff-directory/biography/clore-marius
[pubs]: https://www.niddk.nih.gov/about-niddk/staff-directory/biography/clore-marius/publications
[section]: https://www.niddk.nih.gov/research-funding/at-niddk/labs-branches/laboratory-chemical-physics#section-of-molecular-and-structural-biophysics
[lab]: https://www.niddk.nih.gov/research-funding/at-niddk/labs-branches/laboratory-chemical-physics
[niddk]: https://www.niddk.nih.gov
[apache]: https://httpd.apache.org
