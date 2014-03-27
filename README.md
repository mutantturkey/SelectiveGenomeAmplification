SelectiveGenomeAmplification
============================

PI: http://brisson.bio.upenn.edu/



## Requirements
To use this you'll need:

 - A unix environment
 - kmer_total_count, a kmer counter available here: http://github.com/mutantturkey/dna-utils/
 - bash or compliant shell.
 
 
## Setup

    git clone git@github.com:mutantturkey/SelectiveGenomeAmplification.git
    cd SelectiveGenomeAmplification
    make
    sudo make install

## Usage Examples
Standard use of (SGA) SelectiveGenomeAmplification is easy. it takes two arguments,
the foreground and background


    SelectiveGenomeAmplification PfalciparumGenome.fasta HumanGenome.fasta;
    less PfalciparumGenome_HumanGenome/final_mers

SGA allows for many tunable parameters, which are all explained in the chart
below.  For user customizable variables, they need to be passed in as
environmental variables like so:

    max_mer_distance=5000 max_select=6 min_mer_range=6 max_mer_range=12 \
    SelectiveGenomeAmplification.sh PfalciparumGenome.fasta half.fasta 

SGA also comes with a easy to use user prompt called SelectiveGenomeAmplificationUI.
It allows for a less expereienced user to use
SGA without issue.

### Running individual steps

By default SelectiveGenomeAmplification runs all four steps, but you can
specify the program to run other steps, like in these examples.

    current_run=run_1 SelectiveGenomeAmplification target.fasta bg.fasta score

    current_run=run_1 SelectiveGenomeAmplification target.fasta bg.fasta select score

    current_run=run_1 SelectiveGenomeAmplification target.fasta bg.fasta 3 4 

valid steps are these:

- count (1)
- filter (2)
- select (3)
- score (4)

This function does not try to be smart, so use it wisely
    

### Manually scoring mer combinations

Users can manually score combinations of mers they choose using the
score\_mers.py script.

    score_mers.py -f foreground.fa -b background.fa -c combination file -o output


The combination file should look like this:

    ACGATATAT TACATAGA TATATATAT ACGTACCAT ATATTA
    AAATTATCAGT ATACATA ATATACAT ATATACATA ACATA
		ATATACATA ATCATGATA CCAGATACATAT

each row is combination to be scored.

## Customizable variables

range of mers, min and max 

variable | default | notes
:---- | :---- | ---- | :----
current\_run | Not Enabled | specify the run you want to run steps on
min\_mer\_range | 6  | minimum mer size to use
max\_mer\_range | 12 | maximum mer size to use 
max\_mer\_distance | 5000 | maximum distance between mers in foreground
output\_directory | $foreground\_$background/ | ex. if fg is Bacillus.fasta and  bg is HumanGenome.fasta then folder would be $PWD/Bacillus.fasta\_HumanGenome\_output.fasta/
counts\_directory | $output\_directory/.tmp | directory for counts directory
tmp\_directory | $output\_directory/.tmp | temporary files directory
max\_melting\_temp | 30° | maximum melting temp of mers
min\_melting\_temp | 0° | minimum melting temp of mers
min\_foreground\_binding\_average | 50000 | elminate mers that appear less frequently than the average  (length of foreground / # of occurances)
max\_select | 15 | maximum number of mers to pick
max\_check | 35  | maximum number of mers to select (check the top #)
ignore\_mers | Not Enabled | mers to explicitly ignore, space seperated ex. ignore\_mers="ACAGTA ACCATAA ATATATAT"
foreground | Not Enabled | path of foreground file
background | Not Enabled | path of background file
max\_consecutive\_binding | 4 | The maxium number of consecutive binding nucleotides in homodimer and heterodimers
fg\_weight | 0 | How much extra weight to give higher frequency mers in fg. see "equations" (between 0 and 1)
primer\_weight | 0 | How much extra weight to give to sets with a higher number of priemrs. (between 0 and 1)

## Equations

Here's what we are using to determine our scoring and selectivity

### Selecivity

Our selectivity is what we use to determine what top $max\_check mers are checked later
on in our scoring function. Currently we use this formula:

By default our fg\_weight is zero. This gives no extra weight to more
frequently occuring mers, but can be set higher with the fg\_weight
environmental variable if you wish to do so.

    hit = abundance of primer X (ex. 'ATGTA') in background

    (foreground hit / background hit) * (foreground hit ^ fg_weight)


### Score function
