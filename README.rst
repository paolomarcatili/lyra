BCR and TCR modelling tool
==========================

0. Prerequisites
   In order to install Lyra, you will need the hmmer package installed
   so that the executables hmmalign and hmmsearch can be located in
   any PATH folder. On Debian-like systems, it can be done using apt:

   sudo apt-get install hmmer

   Please refer to the hmmer documentation for additional information.
   It is adviced to install Scwrl4 as well. The Scwrl executable
   should be present in the PATH
   
   
1. Installation

   [optional] It is higly adviced to use python virtual environment:
   cd your_VE_folder
   virtualenv lyra
   source ./lyra/bin/activate
   
   you can install Lyra by using pip on the uncompressed package:

   tar -xzvf lyra.tar.gz
   pip install bcr_models

   It will install all the python required packages. To uninstall, use
   pip uninstall bcr_models

   2. Usage
      The package will install the executable lyra_model. This
      executable allows the modeling of TCRs and BCRs given their
      sequence alone. The basic usage is the following:

      lyra_model sequence.fasta -o model.pdb

      where sequence.fasta is a fasta format file containing the
      sequences of both chains as either nucleotide or amino acid, and
      model.pdb is the output 3d model in pdb format. The model is
      built using automatically selected templates. A structural check
      is performed to detect any interruption or abnormal bond length
      between the model main chain atoms, and the result are printed
      when running in the verbose mode (-v).

      The additional following options are available:

      usage: lyra_model [-h] [-o OUTPUT] [--no-renumber] [-q OUTPUT_NOSCWRL]
                  [-p OUTPUT_PYMOL] [--no-scwrl] [--verbose] [-a] [--debug]
                  [--align] [--scwrl-method {different,all}]
                  seq [hmms [hmms ...]]

      positional arguments:
      seq                   Light and heavy chain sequence in fasta file
      hmms                  HMM profiles to align. Omit to use built-in profiles

      optional arguments:
      -h, --help            show this help message and exit
      -o OUTPUT, --output OUTPUT
                        Output PDB file
	--no-renumber         Renumber the sequences using internal HMM numbering
                        rather that KC numbering scheme
	-q OUTPUT_NOSCWRL, --output-noscwrl OUTPUT_NOSCWRL
                        Output PDB file without running Scwrl
	-p OUTPUT_PYMOL, --output-pymol OUTPUT_PYMOL
                        Output PyMol script
	--no-scwrl            Do not use SCWRL4
	--no-check            Do not perform structural check on the model
	--verbose
	-a, --noalign         Use the alignment provided by the user in the input
	                              file
	--debug               Make an effort to keep things working.
	--align               Only align the sequences
	--scwrl-method {different,all}
                        Scwrl4 method to remodel sidochains: different or all

	If the -a option is used, the input file should contain the aligned
	amino acid sequences according to the internal HMM numbering
	scheme. Be aware that uncorrect alignments might result in
	incomplete models.
	
