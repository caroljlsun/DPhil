#!/bin/bash
COUNTER=0
for OUTPUT in $(seq 5)
do
    let COUNTER++
done
printf "The value of the counter is COUNTER=%d\n" $COUNTER