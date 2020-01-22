# pkin-filter

Simple Script to analyse docking results from mutated versions of protein Kinases  
Script uses slightly modified version of renumber_pdb.py from https://github.com/gawells/pdb_renumber  
Dependencies:  
Prody, Emboss, openbabel (python package, including pybel) and BioPython  
usage: hbondfilter.py [-h] -d DIRECTORY -r REFSEQ_FOLDER -t TARGET -i INPUT

                      [-s] -a APPENDIX [-p]

optional arguments:

  -h, --help            show this help message and exit

  -d DIRECTORY, --directory DIRECTORY

                        Specify folder containing the docking results

  -r REFSEQ_FOLDER, --refseq_folder REFSEQ_FOLDER

                        Specify folder containing the canonical sequences of

                        the Kinases

  -t TARGET, --target TARGET

                        Specify your protein target e.g. EGFR

  -i INPUT, --input INPUT

                        Specify (in lower case) the 3-letter input file format

                        of the docked molecules e.g. sdf

  -s, --sort            use flag if you want to sort data by most passed poses

  -a APPENDIX, --appendix APPENDIX

                        Specify suffix used to differentiate docked molecules

                        from non-docked e.g. _dock.sdf

  -p, --pickled         Run this flag after analysing the data at least once,

                        this flag will load in the previously analysed data

                        without rerunning the analysis
