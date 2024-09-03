# :hatching_chick::dna: Vertebrate age at maturity prediction from genomic data

Analysis pipeline for:

>Budd, A., Yong, S.Y., Heydenrych, M., Mayne, B., Berry, O. and Jarman, S., 2024. Universal prediction of vertebrate species age at maturity.

[Pre-print](https://assets-eu.researchsquare.com/files/rs-4448358/v1/6a612c55-fe0d-46a6-ac3f-b55bf3c22ba7.pdf)

[![DOI](https://zenodo.org/badge/781254240.svg)](https://zenodo.org/doi/10.5281/zenodo.13637779)

---

## :file_folder: Structure

The code is separated into the following three folders:

`A_getGenomicData` - Contains all steps run remotely. Running these is optional (see **Data** section below for details).

`B_exploreData` - Data exploration and partitioning steps. 

`C_model` - All modelling steps and data including figure generation.

`life_pi` - Contains prediction interval steps from https://bitbucket.csiro.au/projects/AI-TOOLBOX_SUKYEE/repos/life_pi/browse and written by Suk Yee Yong (sukyee.yong@csiro.au). Included here for use in figure generation.  

## :chart_with_upwards_trend: Data

All required files are in the repo, with the exception of the raw genomic data which can be downloaded directly from https://www.ncbi.nlm.nih.gov/genome/.

All files and folders should be run in order, however, the code in `A_getGenomicData` was written for and run on a HPC facility that uses a SLURM batch-queue system. This means that many of the slurm scripts (.slurm extension) specify core allocation, run times and memory usage allocation that may need to be adapted for different platforms. To simplify use, all output from the steps in `A_getGenomicData` are included in the repo so that users can run the local steps (in `B_exploreData` and `C_model`) of the analyses *without* having to repeat any of the steps run remotely on a HPC.

Output files begin with the number of the script in which they were generated. **NB: Large files have been compressed and will need to be expanded before use.** 

## :woman_technologist: Author
Alyssa Budd (alyssa.budd@csiro.au).

## :bouquet: Acknowledgements
The scripts have been adapted from a previous model used to predict fish lifespan https://github.com/dr-budd/fish_life. All scripts in `lifepi` were by written by Suk Yee Yong (sukyee.yong@csiro.au).

## :copyright: License
[CSIRO Open Source Software Licence Agreement (variation of the BSD / MIT License)](LICENSE.txt)
