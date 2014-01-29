SelectiveGenomeAmplification
============================

PI: http://brisson.bio.upenn.edu/

To use this you'll need:

 - A unix environment
 - kmer_total_count, a kmer counter available here: http://github.com/mutantturkey/dna-utils/
 - bash or compliant shell.
 
 
Setup:

    git clone git@github.com:mutantturkey/SelectiveGenomeAmplification.gi
    cd SelectiveGenomeAmplification
    make
    sudo make install

Example Usage:

    selectiveGenomeAmplification PfalciparumGenome.fasta HumanGenome.fasta;
    less PfalciparumGenome_HumanGenome/final_mers

For user customizable variables:

    max_mer_distance=5000; max_select=6 min_mer_range=6 max_mer_range=12 \
    SelectiveGenomeAmplification.sh PfalciparumGenome.fasta half.fasta 

## Customizable variables

range of mers, min and max 

C | variable | default | notes
:---- | :---- | :---- | ---- | :----
Y | min_mer_range | 6  | minimum mer size to use
Y | max_mer_range | 10 | maximum mer size to use 
Y | max_mer_distance | 5000 | maximum distance between mers in foreground
N | output_directory | $PWD/$foreground_$background/ | ex. if fg is Bacillus.fasta and  bg is HumanGenome.fasta then folder would be $PWD/Bacillus.fasta_HumanGenome_output.fasta/
Y | counts_directory | $output_directory/.tmp | directory for counts directory
Y | tmp_directory=$output_directory/.tmp | temporary files directory
Y | max_melting_temp | 30° | maximum melting temp of mers
Y | min_melting_temp | 0° | minimum melting temp of mers
Y | min_mer_count | Not Enabled (0) | only select mers that occur more frequently than this number
Y | max\_select | maximum number of mers to pick
Y | ignore\_mers | Not Enabled | mers to explicitly ignore, space seperated ex. ignore\_mers="ACAGTA ACCATAA ATATATAT"
Y | foreground | Not Enabled | path of foreground file
Y | background | Not Enabled | path of background file
Y | max\_consecutive\_binding | 4 | The maxium number of consecutive binding nucleotides in homodimer and heterodimers

