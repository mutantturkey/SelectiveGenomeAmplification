#!/usr/bin/env bash


function clean_exit {

	kill -s 9 $pid
	exit
}
trap clean_exit SIGHUP SIGINT


echo "score_mers.py $1 $2 $3 $4 &"
score_mers.py $1 $2 $3 $4 
pid=$!

wait $pid
