SelectiveWholeGenomeAmplification
============================

SWGA is a tool for choosing primers for the selective amplification of a
target genome from a sample containing a mixture of target and contaminating
DNA (i.e. pathogen genome from infected host blood) [cite relevant paper]. It
does so by identifying short, recurring motifs in a target sequence file and
scoring sets of motifs based on selectivity for and even distribution in the
target sequence against a background sequence file.

PI: http://brisson.bio.upenn.edu/

## Table of Contents

* [Requirements](#requirements)
* [Setup](#setup)
* [Example Usage](#example-usage)
  * [SWGA User Interface](#sga-user-interface)
  * [Setting Tunable Parameters](#setting-tunable-parameters)
  * [Running Individual Steps](#running-individual-steps)
  * [Manually Scoring Specific Mer Combinations From List ](#manually-scoring-specific-mer-combinations-from-list)
  * [Manually Score All Combinations From List](#manually-score-all-combinations-from-list)
  * [Manually Rescore All Combinations From Previously Scored File](#manually-rescore-all-combinations-from-previously-scored-file)
* [Table of Tunable Parameters](#tunable-parameters)
* [Equations](#equations)
  * [Mer Selectivity](#mer-selectivity)
  * [Scoring Combinations](#score-combinations)
    * [Default Scoring Function](#default-scoring-function)
    * [Custom Scoring Function](#custom-scoring-function)
* [Filters](#filters)
* [Output](#output)
  * [Select Mers](#select_merspy-output)
  * [Score Mers](#score_merspy-output)
  
## Requirements
To use this you'll need:

 - A Unix environment
 - [dna-utils](http://github.com/mutantturkey/dna-utils/)
 - bash or compliant shell.
 - python 2.7.x
 
## Setup

    git clone git@github.com:mutantturkey/SelectiveWholeGenomeAmplification.git
    cd SelectiveWholeGenomeAmplification
    make
    sudo make install

## Example Usage
Standard use of (SGA) SelectiveWholeGenomeAmplification is easy. it takes two arguments,
the foreground and background


    SelectiveWholeGenomeAmplification PfalciparumGenome.fasta HumanGenome.fasta;
    less PfalciparumGenome_HumanGenome/final_mers

### SWGA User Interface
SWGA also comes with a easy to use user prompt called SelectiveWholeGenomeAmplificationUI.
It allows for a less experienced user to use SWGA without issue. to run this
all you need to do is run SelectiveGenomeAmiplifcationUI and you'll see a
series of prompts asking the user about tunables like below

    Where would you like to temporary files to be stored? (Default=$output_directory/.tmp): 
    Where would you like to count files to be stored? (Default=$output_directory/.tmp): 
    maximum mer size you would like to pick? (Default=12): 10
    minimum mer size you would like to pick? (Default=6): 7
    eliminate mers that appear less frequently on average than this number ? (Default=50000): 25000
    .....
    Input the path to your foreground file:target.fa  
    Input the path to your background file:humangenome.fa 
    Would you like to output your inserted variables to a string you can later paste? (Y/N/Default=y): n
    Run SelectiveWholeGenomeAmplification? (Y/N/Default=y): y

### Setting Tunable Parameters

SGA allows for many tunable parameters, which are all explained in the chart
below.  For user customizable variables, they need to be passed in as
environmental variables like so:

    max_mer_distance=5000 max_select=6 min_mer_range=6 max_mer_range=12 \
    SelectiveWholeGenomeAmplification.sh PfalciparumGenome.fasta half.fasta 


### Running individual steps

By default SelectiveWholeGenomeAmplification runs all four steps, but you can
specify the program to run other steps, like in these examples.

    current_run=run_1 SelectiveWholeGenomeAmplification target.fasta bg.fasta score

    current_run=run_1 SelectiveWholeGenomeAmplification target.fasta bg.fasta select score

    current_run=run_1 SelectiveWholeGenomeAmplification target.fasta bg.fasta 3 4 

valid steps are these:

- count (1)
- filter (2)
- select (3)
- score (4)

This function does not try to be smart, so use it wisely.

### Manually scoring specific mer combinations from list

Users can manually score combinations of mers they choose using the
score\_mers.py script.

    score_mers.py -f foreground.fa -b background.fa -c combination file -o output


The combination file should look like this:

    ACGATATAT TACATAGA TATATATAT ACGTACCAT ATATTA
    AAATTATCAGT ATACATA ATATACAT ATATACATA ACATA
    ATATACATA ATCATGATA CCAGATACATAT

each row is combination to be scored.


### Manually score all combinations from list

Users can manually score all  combinations of mers they choose using the
score\_mers.py script.

    score_mers.py -f foreground.fa -b background.fa -m mer file -o output


The mer file should look like this:

    ATATAT
    TACATA
    TACATAGCA
    TATAGAATAC
    CGTAGATA
    TAGAAT

each row is a separate mer. do not put multiple mers on one line.

### Manually rescore all combinations from previously scored file

Users can manually rescore all combinations of mers they previously used in the
score\_mers.py script. This allows users to test different score functions easily
with the same combinations.

An example would be this:

    score_func=nb_primers**2 score_mers.py -f fg.fa -b bg.fa -r fg_bg/run_1/all-scores -o primers_squared_scores


## Tunable Parameters

variable | default | notes
:---- | :---- | ---- | :----
current\_run | Not Enabled | specify the run you want to run steps on
min\_mer\_range | 6  | minimum mer size to use
max\_mer\_range | 12 | maximum mer size to use 
max\_mer\_distance | 5000 | maximum distance between mers in foreground
min\_melting\_temp | 0° | minimum melting temp of mers
max\_melting\_temp | 30° | maximum melting temp of mers
min\_foreground\_binding\_average | 50000 | eliminate mers that appear less frequently than the average  (length of foreground / # of occurrances)
ignore\_mers | Not Enabled | mers to explicitly ignore, space separated ex. ignore\_mers="ACAGTA ACCATAA ATATATAT"
ignore\_all\_mers\_from\_files | Not Enabled | ignore any mers found in these files. space separated.
output\_directory | $foreground\_$background/ | ex. if fg is Bacillus.fasta and  bg is HumanGenome.fasta then folder would be $PWD/Bacillus.fasta\_HumanGenome\_output.fasta/
counts\_directory | $output\_directory/.tmp | directory for counts directory
tmp\_directory | $output\_directory/.tmp | temporary files directory
max\_select | 15 | maximum number of mers to pick
max\_check | 35  | maximum number of mers to select (check the top #)
foreground | Not Enabled | path of foreground file
background | Not Enabled | path of background file
max\_consecutive\_binding | 4 | The maximum number of consecutive binding nucleotides in homodimer and heterodimers
fg\_weight | 0 | How much extra weight to give higher frequency mers in fg. see "equations" (between 0 and 1)
primer\_weight | 0 | How much extra weight to give to sets with a higher number of primers. (between 0 and 1)
output\_top\_nb | 10000 | How many scores do you want to output in your sorted output file?
score\_func | Not Enabled | see the [custom scoring](#custom-scoring-function) section
sort\_by | min | How do you want to rank top-scores? min means smaller is better, max is larger. 'min' or 'max'

## Equations

Here's what we are using to determine our scoring and selectivity

### Mer Selectivity

Our selectivity is what we use to determine what top $max\_check mers are checked later
on in our scoring function. Currently we use this formula:

By default our fg\_weight is zero. This gives no extra weight to more
frequently occurring mers, but can be set higher with the fg\_weight
environmental variable if you wish to do so.

    hit = abundance of primer X (ex. 'ATGTA') in background

    (foreground hit / background hit) * (foreground hit ^ fg_weight)


### Scoring combinations 

All variables used in our scoring function are described here:

    fg_pts = an array of all the points of each mer in the combination, and sequence ends
    fg_mean_dist = mean distance between each point in fg_pts
    fg_stddev = standard deviation of distance between each point in fg_pts

    nb_primers = number of primers in a combination
    primer_weight = extra weight for sets with higher primers

    bg_ratio = length of background / number of times primer was in background

#### Default scoring function

The default scoring function is this:

    mer_score = (nb_primers**primer_weight) * (fg_mean_dist * fg_std_dist) / bg_ratio

#### Custom scoring function

We support custom scoring via python's exec methods. This means that you can
destroy your system, blow up the universe, implode your hard drive, all within
the confines of this exec. That means don't do anything crazy. Stick to basic arithmetic.
 
This is a security hole.

you can specify it like any other parameter like so:

    # the default function
    score_func="(nb_primers**primer_weight) * (fg_mean_dist * fg_std_dist) / bg_ratio"

You need to use **valid** python code. 

## Filters

There are several filters that our mers go through, to eliminate ones that won't fit our needs. They are all configurable via the tunable parameters. If you look in a output directory, you'll see a folder called "passes-filter". This contains a file for each of the different steps in the pipeline, and the contents of each file is what 'passes' that filter.

For example, if you ignored the mer 'AAAAA', then in passes-filter/1-$foreground-ignore-mers there would be no line containing that.

The filter system works like a big pipe, whatever gets filtered out won't make it to the next step. the order is like this


    All mers -> ignore_mers -> ignore_all_mers -> average_binding -> non_melting -> consecutive_binding
## Output

The file structure outputted by default is this:

    $foreground_$background
    └── run_1 # current_run
        ├── passes-filter # filter folder for filtering steps
        │   ├── 1-$foreground-ignore-mers
        │   ├── 2-$foreground-ignore-all-mers
        │   ├── 3-$foreground-average-binding
        │   ├── 4-$foreground-non-melting
        │   └── 5-$foreground-consecutive-binding
        ├── $foreground-filtered-counts # final filtered mers used for select_mers.py
        ├── parameters # parameters used in the run
        ├── selected-mers # final filtered mers used for select_mers.py
        ├── selected-mers # final filtered mers used for select_mers.py
        ├── all-scores    # file outputted by score_mers.py (all the scores generated)
        └── top-scores    # the sorted top $output_top_nb scores from all-scores

### select\_mers.py output

Select mers outputs a tab delimited file, with 4 columns: mer, foreground count,
background count, and the mer selectivity value. (higher is better)

    CTAACTTAGGTC  1572  155  10.14194
    CTAACATAGGTC  1479  132  11.20455
    GACCTATGTTAG  1479  132  11.20455


### score\_mers.py output

score mers outputs a tab delimited file with 6 columns:

    nb_primers  Combination  Score  FG_mean_dist  FG_stdev_dist  BG_ratio
