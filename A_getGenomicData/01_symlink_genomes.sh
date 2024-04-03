#!/bin/bash

## SET UP ----

## create variables for directories and files

WD=$PWD
DD=$WD/dataFiles
GD=$WD/../chordGenomeDb/chordGenome
OD=$WD/SymVertGenomes

## make the directories if they don't already exist

if [ ! -d "$OD" ]
 then
	mkdir -p "$OD"
fi

## START SCRIPT ----

## NB this script symlinks genomes from an existing database to the working directory
## if you don't already have the genomes, you will need to download them from NCBI
## using the accession numbers provided in the supplementary material

## read the text file line by line
while IFS=$'\t' read -r -a line; do

    ## get the file names
    GENOME="${line[0]}"
    CHECKSUMS="${line[1]}"
    SPECIES="${line[2]}"

    ## create the symlinks in the destination folder if required
	if [ -e $OD/$GENOME ]; then
		echo $SPECIES" genome ("$GENOME") is already symlinked"
	else
		echo $SPECIES" genome ("$GENOME") and checksums ** symlinking now **"
		ln -s "$GD/$GENOME" "$OD/$GENOME"
		ln -s "$GD/$CHECKSUMS" "$OD/$CHECKSUMS"
	fi

done < $DD/00.03_vert_agemat_genome_symlist.txt > 01_symlink_genomes.log

## END SCRIPT
