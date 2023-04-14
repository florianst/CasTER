# CasTER
CasTER (Cas9 Thermodynamic Energy Repair prediction) is a tool to predict the outcome of cell-internal repair machinery after Cas9-induced DNA double strand breaks. It uses a tree-based model based on thermodynamic features generated using [MOSAICS](https://www.cs.ox.ac.uk/mosaics/Downloads.php). The automatic generation of PDB files from DNA sequences is done using [fd_helix.c](https://casegroup.rutgers.edu/).

CasTER is developed by [Florian Störtz](https://users.ox.ac.uk/~lina3430/) at the group of Dr. Minary, University of Oxford.

## Usage
Running CasTER requires two steps:
```
python load_data.py
```
Positional arguments are:
- oligo file: path to input file of tab separated oligos, containing the following columns: ID, Guide Sequence, Target Sequence, Scaffold, Subset, Comments, PAM Index, Strand. Only ID, Target Sequence, PAM Index and Strand need to have meaningful content. For an example, see ```oligos.txt```
- number of CPU workers (default: 48 - choose this proportional to the number of CPUs available)
- width in base pairs around the cut site which is rendered and used by MOSAICS (default: 21)
- indel file: path to a file of measured indels, each oligo introduced by ```@@@ID```, then following an enumeration of indels with tab-separated columns: indel signature, read count, repaired sequence. For an example, see ```indels.txt```. This file does not need to be supplied, but if it is supplied then CasTER calculates benchmarks (KL divergence, Pearson r, Spearman r) based on the difference between true and predicted indel distributions.

This generates the ```mosaics_data.pickle``` file containing the MOSAICS energies which are needed for prediction in the following step.

To run a prediction, call
```
python caster.py
```
Positional arguments are:
- sum_level: ```0```: predict indel length, ```1```: predict indel position, ```None```: predict indel length and position
- feature_set: ```0```: use sequence and energy features, ```1```: use energy features only, ```2```: use sequence features only
- mosaics_data_path: path to the file containing the MOSAICS energies (default: ```mosaics_data.pickle``` )

### Example
Running ```python caster.py``` gives you an example output using precalculated MOSAICS energy features on 

## Installation
Download [MOSAICS](https://www.cs.ox.ac.uk/mosaics/Downloads.php) and install using the [instructions](https://www.cs.ox.ac.uk/mosaics/Downloads.php#standaloneinstall) given there. Move ```mosaics.x``` to the CasTER directory.

Adjust the paths to the MOSAICS library folder in ```input```

### Requirements
See ```requirements.txt```

