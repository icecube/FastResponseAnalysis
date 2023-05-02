#!/bin/bash

while true
do
    echo 
    echo $(date -u) Starting to update plots
    python Data_Display.py
    echo $(date -u) Plots Updated OK
	sleep 3600
done