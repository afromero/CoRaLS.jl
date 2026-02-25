#!/usr/bin/env bash
#SBATCH -A PAS0654
#SBATCH --job-name=accpt_array
#SBATCH --mail-type=BEGIN,END,FAIL ## This is so you get an email of when the job starts and finishes
#SBATCH --output=out/corals_%A_50km_%am.out
#SBATCH --error=err/corals_%A_50km_%am.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem=4G

## Example: sbatch --array=1-20 --export=ALL,ENERGY=1,ICE=5 array_altitude.sh
## Example: sbatch --array=1-20 --export=ALL,ALT=1,ENERGY=1,ICE=5,ANT=4,TRIG=4,ANG=-90,FREQ1=300,TEXP=7,VAR=ALT array_altitude.sh

#cd ~/../../../fs/scratch/PAS2277/linton93/CoRaLS_MC/

cd ~/CoRaLS/main/ ## If using the main branch
cd ~/CoRaLS/dev/ ## If using a development branch

## Set the variables we read in

#ALT=$1					# altitude in km
#ENERGY=$2			# Energy multiplier (set to 1)
#ICE=$3					# ice depth in m
#ANT=$4					# number of antennas
#TRIG=$5				# number of triggers required
#ANG=$6					# pointing angle of antennas in degrees
#FREQ1=$7				# minimum frequency in band in MHz
#TEXP=$8				# exponent on 10 for ntrials


# Slurm‐provided vars
TASK="$SLURM_ARRAY_TASK_ID"

## Whichever of the above variables you want to vary with task
##  set equal to $TASK here:
eval ${VAR}=${TASK}
echo $VAR
echo ${ALT}

#if [[ $VAR == "ALT" ]]
#then
#				echo "alt=$((5 * ${ALT})) km   energyMult=${ENERGY}   ice=${ICE} m   ant=${ANT}   trig=${TRIG}   angle=${ANG} deg  freqMin=${FREQ1} MHz   TEXP=${TEXP}"
#elif [[ $VAR == "ICE" ]]
#then
#
#				echo "ice=${ICE} m   energyMult=${ENERGY}   alt=$((5 * ${ALT})) km   ant=${ANT}   trig=${TRIG}   angle=${ANG} deg  freqMin=${FREQ1} MHz   TEXP=${TEXP}"
#fi
## call Julia with (altitude, ice_depth, bin_start, bin_end)
#julia CoRaLS.jl/slurm/generic_acceptance.jl "$ALT" "$ENERGY" "$ICE" "$ANT" "$TRIG" "$ANG" "$FREQ1" "$TEXP" 


if [[ $VAR == "ALT" ]]
then
        echo "Varying altitude"
        echo "alt=$((5 * ${ALT})) km   energyMult=${ENERGY}   ice=${ICE} m   ant=${ANT}   trig=${TRIG}   angle=${ANG} deg  freqMin=${FREQ1} MHz   TEXP=${TEXP}"
elif [[ $VAR == "ICE" ]]
then
        echo "Varying ice depth"
        echo "alt=$((5 * ${ALT})) km   energyMult=${ENERGY}   ice=${ICE} m   ant=${ANT}   trig=${TRIG}   angle=${ANG} deg  freqMin=${FREQ1} MHz   TEXP=${TEXP}"
elif [[ $VAR == "ANG" ]]
then
        echo "Varying antenna angle"
        echo "alt=$((5 * ${ALT})) km   energyMult=${ENERGY}   ice=${ICE} m   ant=${ANT}   trig=${TRIG}   angle=${ANG} deg  freqMin=${FREQ1} MHz   TEXP=${TEXP}"
fi

# call Julia with (altitude, ice_depth, bin_start, bin_end)
julia CoRaLS.jl/slurm/generic_acceptance.jl "$ALT" "$ENERGY" "$ICE" "$ANT" "$TRIG" "$ANG" "$FREQ1" "$TEXP"



#echo "Done!"
