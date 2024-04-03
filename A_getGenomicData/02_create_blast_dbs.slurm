#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=72:00:00
#SBATCH --output=%x.%jx.log

## SET UP ----

## create variables for directories and files
WD=$PWD
DD=$WD/SymVertGenomes
OD=$WD/BlastDatabases

## make the output directory if it doesn't already exist
if [ ! -d "$OD" ]
 then
    mkdir -p "$OD"
fi

## load module(s)
module load blast+/2.12.0

## RUN ----

## print date time start
echo -n "start: " && date

for FILE in $DD/*.fna.gz
do
  NAME=$(basename $FILE)
  ID=${NAME%.fna.gz}
  DBNAME=$(echo $ID'_db')
  DBFILE=$(echo $DBNAME'.ndb')
  if [ -e $OD/$DBFILE ]; then
    echo $DBNAME" already created"
  else
    echo "*** creating "$DBNAME
    ## decompress the file and pass the output to makeblastdb
    gunzip -c "$FILE" | makeblastdb -dbtype nucl -in - -out "$OD/$DBNAME" -title "$ID"
  fi
done

## print date time end
echo -n "end: " && date

## END SCRIPT
