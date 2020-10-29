#!/bin/bash
K=5
Z=1
Y=3
java LearnTopicModel -model flda -input ~/flda_res/chinese_1000_place.txt -K $K -Z $Z -Y $Y -iters 5000 -sigmaW 0.1
java LearnTopicModel -model flda -input ~/flda_res/chinese_1000_location.txt -K $K -Z $Z -Y $Y -iters 5000 -sigmaW 0.1
java LearnTopicModel -model flda -input ~/flda_res/chinese_1000_service.txt -K $K -Z $Z -Y $Y -iters 5000 -sigmaW 0.1
python2 topwords_flda.py $K $Z $Y ~/flda_res/chinese_1000_place.txt 100 > ~/flda_res/chinese_place_output.txt
python2 topwords_flda.py $K $Z $Y ~/flda_res/chinese_1000_location.txt 100 > ~/flda_res/chinese_location_output.txt
python2 topwords_flda.py $K $Z $Y ~/flda_res/chinese_1000_service.txt 100 > ~/flda_res/chinese_service_output.txt
