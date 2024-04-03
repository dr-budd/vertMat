#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32GB
#SBATCH --time=08:00:00
#SBATCH --output=%x.%j.log

## SET UP ----

## create variables for directories and files

WD=$PWD
DD=$WD/BlastDatabases

MODH=$WD/QueryResultsMammalsHomo
FOD=$WD/QueryResultsFishDanio
ROD=$WD/QueryResultsReptilesGallus

HUMAN_PROMOTERS=$WD/promotersCore/epd_homo_sapiens_most_rep_core.fa
ZFISH_PROMOTERS=$WD/promotersCore/epd_danio_rerio_most_rep_core.fa
CHICK_PROMOTERS=$WD/promotersCore/epd_gallus_gallus_most_rep_core.fa

## remove old and create new directories
rm -r $MODH $MODM $FOD $ROD
mkdir $MODH $MODM $FOD $ROD

## load module(s)
module load blast+/2.12.0

## RUN ----

## print date time start
echo -n "start: " && date

## read the text file line by line
while IFS=$'\t' read -r -a line; do

    ## get the file names
    GENOME="${line[0]}"
    CHECKSUMS="${line[1]}"
    SPECIES="${line[2]}"
    GROUP="${line[4]}"
    ID=$(basename $GENOME .fna.gz)

        DBNAME=$(echo $ID'_db')
        THNAME=$(echo $ID'_top_hits')

        # echo $DBNAME $THNAME

        if [ $GROUP == "Reptiles" ]; then
                echo "querying chicken promoters for "$DBNAME" for the following reptile (inc. birds): "$SPECIES
                blastn -query $CHICK_PROMOTERS -db $DD/$DBNAME -out $ROD/$THNAME -outfmt "6 qseqid sseqid slen qstart qend sstart send pident nident sseq" -perc_identity 70  -culling_limit 1
        else
                if [ $GROUP == "Mammals" ]; then
                        echo "querying human promoters for "$DBNAME" for the following mammal: "$SPECIES
                        blastn -query $HUMAN_PROMOTERS -db $DD/$DBNAME -out $MODH/$THNAME -outfmt "6 qseqid sseqid slen qstart qend sstart send pident nident sseq" -perc_identity 70  -culling_limit 1
                else
                        if [ $GROUP == "Fish" ]; then
                                echo "querying zebrafish promoters "$DBNAME" for the following fish: "$SPECIES
                                blastn -query $ZFISH_PROMOTERS -db $DD/$DBNAME -out $FOD/$THNAME -outfmt "6 qseqid sseqid slen qstart qend sstart send pident nident sseq" -perc_identity 70  -culling_limit 1
                        else echo "*** $SPECIES is not a mammal, bird/reptile or fish. No query made."

                        fi
                fi
        fi

done < $WD/dataFiles/00.03_vert_agemat_genome_symlist.txt

## print date time end
echo -n "end: " && date

## END SCRIPT