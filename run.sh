#!/bin/bash
K=$1
Z=$2
Y=$3
inputFile="/home/ubuntu/data/total/chinese_total.txt"

java LearnTopicModel -model flda -input $inputFile -K $K -Z $Z -Y $Y -iters 5000 -sigmaW 0.1
python2 python2 topwords_flda.py $inputFile 100 > ~/data/total/chinese_total_output.txt
