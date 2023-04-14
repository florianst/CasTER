#!/bin/bash

trim() {
    local var="$*"
    # remove leading whitespace characters
    var="${var#"${var%%[![:space:]]*}"}"
    # remove trailing whitespace characters
    var="${var%"${var##*[![:space:]]}"}"   
    echo -n "$var"
}

FILENAME="run.pdb"
TEMPFILE="temp.pdb"

if [ $# -eq 3 ]; then
    SEQUENCE=$1
    BPPT=$2
    REMBASESA=0
    REMBASESB=0
    INPUTFILE=$3
elif [ $# -eq 5 ]; then
	SEQUENCE=$1
    BPPT=$2
    REMBASESA=0
    REMBASESB=0
    INPUTFILE=$3
    FILENAME=$4
    TEMPFILE=$5
elif [ $# -eq 7 ]; then
	SEQUENCE=$1
    BPPT=$2
    REMBASESA=$3
    # how many bases to remove from beginning of strand A
    REMBASESB=$4
    # how many bases to remove from end of strand B
    INPUTFILE=$5
    FILENAME=$6
    TEMPFILE=$7
else
    SEQUENCE="AGGCTTAGCGCGGTGTATTGATACCCGCGACTACGGTCCGGGGGGGCTTTGCTAAGCGGT"
    BPPT=10.5
    REMBASES=0
    INPUTFILE="input"
fi

sequence_len=${#SEQUENCE}
sequence_len=$((sequence_len+4))

if [[ $REMBASESA -lt 0 ]]; then
    REMBASESA=0
fi

if [[ $REMBASESB -lt 0 ]]; then
    REMBASESB=0
fi

# generate pdb file from sequence and number of base pairs per turn
# pad the sequence with two bp on each end so we make sure we have end groups on each side after deleting them again
"./fd_helix" "abdna" "AA${SEQUENCE}AA" "$BPPT" > "$FILENAME"

SCRIPTPATH="~/mosaics/MOSAICS_scripts_database/pdb_handlers"
# remove hydrogen atoms (not strictly necessary for 3DDART pdb files)
#(exec "${SCRIPTPATH}/remove_H.pl" "$FILENAME" > "/dev/null")

# add hydrogens using pymol - don't show pymol output
pymol -c "$FILENAME" -d "h_add; save $FILENAME; quit" > "/dev/null"

# find where residue 3 of strand A begins
numline=1 # lines seem to be one-based in sed
residue_firstA=1
while read line; do
    # remove all spaces - desired residue number could be one or two digits!
	residueNum=`echo ${line:23:3} | sed -e "s/^[[:space:]]*//"`
    if [[ $residueNum == $((REMBASESA+3)) ]]; then
    	residue_firstA=$((numline-1))
    	break
    fi
    numline=$((numline+1))
done <$FILENAME

# find where residue sequence_len-1 of strand A begins
numline=1 # lines seem to be one-based in sed
residue_lastA=1
while read line; do
	residueNum=`echo ${line:23:3} | sed -e "s/^[[:space:]]*//"`
    if [[ $residueNum == $((sequence_len-1)) ]]; then
        residue_lastA=$numline
        break
    fi
    numline=$((numline+1))
done <$FILENAME

# find where residue 3 of strand B begins
numline=1 # lines seem to be one-based in sed
residue_firstB=1
while read line; do
	residueNum=`echo ${line:23:3} | sed -e "s/^[[:space:]]*//"`
    if [[ $residueNum == $((sequence_len+3)) ]]; then
    	residue_firstB=$((numline-1))
    	break
    fi
    numline=$((numline+1))
done <$FILENAME

# find where residue sequence_len-1 of strand B begins
numline=1 # lines seem to be one-based in sed
residue_lastB=1
while read line; do
	residueNum=`echo ${line:23:3} | sed -e "s/^[[:space:]]*//"`
    if [[ $residueNum == $(((sequence_len*2)-REMBASESB-1)) ]]; then
    	residue_lastB=$((numline-1))
    	break
    fi
    numline=$((numline+1))
done <$FILENAME

# delete residues which stem from added two bases on each side (for both strands)
sed -i -e 1,${residue_firstA}d "$FILENAME"
residue_lastA=$((residue_lastA-residue_firstA)) # correct line numbers for deletion
residue_firstB=$((residue_firstB-residue_firstA))
residue_lastB=$((residue_lastB-residue_firstA))

sed -i -e $residue_lastA,${residue_firstB}d "$FILENAME"   
residue_lastB=$((residue_lastB-(residue_firstB-residue_lastA))) # correct line number for deletion

sed -i "${residue_lastB},\$d" "$FILENAME"    

# rename base pairs
(exec "${SCRIPTPATH}/rename_DNA_res.pl" "$FILENAME" > "$TEMPFILE")
mv "$TEMPFILE" "$FILENAME"

# correct hydrogen naming
(exec "${SCRIPTPATH}/rename_DNA_H.pl" "$FILENAME" > "$TEMPFILE")
mv "$TEMPFILE" "$FILENAME"

# convert to PAS pdb format
(exec "${SCRIPTPATH}/pdb_to_PAS_pdb.pl" "$FILENAME" > "$TEMPFILE")
mv "$TEMPFILE" "$FILENAME"

# finally, run MOSAICS with correct pos_init_file, eighess_file etc
# TODO: run short cg optimisation before running samc - use out.pdb as input file for samc
("./mosaics.x" "$INPUTFILE")
