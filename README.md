01/05/25 - Simulating Kaon Production in Neutrino Interactions -

This repository contains all of the files from our semester 7 project work and all of the preliminary/support files
for our semester 8 project work. The actual main implementation work can be found at: https://github.com/mmillns2/nuwro

This repo contains 5 main folders; 
- "simulation" and "testing", which contain our very early work from semester 7 and all of the .cc files for the Matrix
  Elements (MEs), as well as testing scripts for these MEs;
- "BASH" contains all of the shell scripts used to run nuwro overnight and generate histograms from the simulated data,
  the "params-files"  folder contains all of the nuwro params files from which the data was generated.
- "mathematica" contains twelve mathematica notebooks which were used for evaluating the MEs from the given hadronic 
  currents. For Sigma and Lambda crossing terms the notebooks labelled 'new' are the up to date ones, the others are old
	testing notebooks.
- "Sem-8-ROOT" contains three sub folders; "all-testing" which contains all our scripts and histograms from our ongoing
  implementations into nuwro throughout semester 8; "scripts" contains all of the final scripts used to generate our
	histograms; "final-plots" contains all of our histograms generated from the simulation data provided by nuwro.
