#!/usr/bin/env bash


function clean_exit {
	echo -e "\nInterrupt caught. Exiting cleanly..."
	echo "killing $pid"
	kill -11 $pid
	sleep .5
	exit
}
trap clean_exit SIGINT

score_mers.py $1 $2 $3 $4 & </dev/null
pid=$!

wait $pid
