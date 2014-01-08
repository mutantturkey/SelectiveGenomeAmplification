#!/bin/bash
# range of mers, min and max 
: ${min_mer_range=6}
: ${max_mer_range=10}
# directory to store our counts and sorted counts
: ${counts_directory=$PWD/counts}
# temp directory 
: ${tmp_directory=$PWD/tmp}
# maximum kmer melting point
: ${max_melting_temp=20}
# minimum mer count
# ! you can supply either a percentage like this .5
# ! or you can supply a raw number (100)
: ${min_mer_count=1000}
# maximum mers to pick
: ${max_select=15}

export min_mer_range
export max_mer_range
export max_select
export min_mer_count


if [ ! -d $counts_directory ]; then
	mkdir $counts_directory
fi

if [ ! -d $tmp_directory ]; then
	mkdir $tmp_directory
fi

foreground=$1
background=$2

if [[ ! -f $foreground ]]; then
	echo "Could not open $foreground."
	exit 1
fi

if [[ ! -f $background ]]; then
	echo "Could not open $background."
	exit 1
fi

for fasta_file in $foreground $background; do

	counts=$counts_directory/$(basename $fasta_file)
	tmp=$tmp_directory/$(basename $fasta_file)

	echo pre-processing $fasta_file

	# check if our preprocessed file exists
	if [[ ! -f $tmp ]]; then
		echo "> pre processed $fasta_file" >> $tmp
		cat $fasta_file | grep -v "^>" | tr -d '\n' >> $tmp
	fi

	# run counts if they haven't been created 
	rm $counts-counts
	for mer in `seq $min_mer_range $max_mer_range`;	do 
		if [ ! -e $counts-counts-$mer ]; then
			echo checking $mer mers for $fasta_file
			kmer_total_count -i $tmp -k $mer -l -n >> $counts-counts-$mer
		else 
			echo "$mer mers already done for $fasta_file"
		fi
		
		cat $counts-counts-$mer >> $counts-counts
	
	done
done


fg_counts=$counts_directory/$(basename $foreground)-counts
bg_counts=$counts_directory/$(basename $background)-counts

fg_tmp=$tmp_directory/$(basename $foreground)
bg_tmp=$tmp_directory/$(basename $background)

echo "checking if mers are below melting temperature in the foreground"

rm $fg_counts-fg-non-melting

python below_melting_temperature.py $max_melting_temp < $fg_counts > $fg_counts-fg-non-melting

python ./select_mers.py $fg_counts-fg-non-melting $fg_tmp $bg_counts $bg_tmp # > $(basename $foreground)_$(basename $background)_final_mers
