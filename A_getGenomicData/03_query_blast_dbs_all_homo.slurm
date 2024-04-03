#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32GB
#SBATCH --time=72:00:00
#SBATCH --output=%x.%j.log

## SET UP ----

## create variables for directories and files
WD=$PWD
DD=$WD/BlastDatabases
GD=$WD/SymVertGenomes
OD=$WD/QueryResultsAllHomo
PROMOTERS=$WD/promotersCore/epd_homo_sapiens_most_rep_core.fa       ## path to promoters

## remove old and make new dir
rm -r $OD
mkdir $OD

## load module(s)
module load blast+/2.12.0

## RUN JOB ----

## print date time start
echo -n "start: " && date

for FILE in $GD/*.fna.gz
do
  NAME=$(basename $FILE)
  ID=${NAME%.fna.gz}
  DBNAME=$(echo $ID'_db')
  THNAME=$(echo $ID'_top_hits')
  if [ -e $OD/$THNAME ]; then
    echo $DBNAME" already queried"
  else
    echo "*** querying "$DBNAME
    blastn -query $PROMOTERS -db $DD/$DBNAME -out $OD/$THNAME -outfmt "6 qseqid sseqid slen qstart qend sstart send pident nident sseq" -perc_identity 70  -culling_limit 1
  fi
done

## print date time end
echo -n "end: " && date

## END SCRIPT
