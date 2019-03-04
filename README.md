# Robotics

This repository contains the localisation code for the kidnapped robot problem. BotSimLib0.37 is created by Austin Gregg-Smith ag7751@bristol.ac.uk and serves as the basis for the code. It contains the main bot simulator as well as a few learning examples. An unworked sample of the library is included. 

The main work is located inside myBotSim folder. localise_trial.m was a longer code made that combined the map and localisation code to make it easier to debug. The final code is localise.m. 

The main method was using a bayesian particle filter to determine which particles are closest to the real position of the lost robot. Path planning was made by using k-nn clustering to make a simplified map and djikstras algorithm. 
