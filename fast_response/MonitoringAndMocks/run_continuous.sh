#!/bin/bash

while true
do
    echo 
    echo $(date -u) Starting to update plots
    python Data_Display.py
    echo $(date -u) Plots Updated OK

    current=$(date -d $(date -d "today" +%H:%M) +%s)
    target=$(date -d $(date -d "$(date -d "today + 1 hour" +%H):30" +%H:%M) +%s)
    sleep $(( $target-$current ))
done