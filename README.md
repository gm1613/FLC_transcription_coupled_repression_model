# FLC_transcription_coupled_repression_model <br>
Code for simulations and analysis used in the manuscript "Proximal termination generates a transcriptional state that determines the rate of establishment of Polycomb silencing "<br>
The file **20230117_autopathway_flc.cpp** contains the source code for simulations of the stochastic model using the Gillespie SSA, written in C++.
The Gillespie SSA implementation uses a pseudo-random number generator from the GNU Scientific Library. The library needs to be linked during compilation.<br>
**The source code was compiled on a system running macOS 12.5, using the clang(LLVM) compiler using the following terminal command:**<br>
clang++ -Wall -pedantic 20230117_autopathway_flc.cpp -o 20230117_autopathway_flc  -lgsl -lgslcblas<br>
<br>
**The resulting output file "20230117_autopathway_flc" can be executed using the following terminal command, including six user specified inputs separated by spaces:**<br>
./20230117_autopathway_flc <"Number of simulated trajectories">  <"Duration of individual simulation(number of cell cycles)"> <"FCA parameter"> <"FLD mediated demethylation probability per histone per proximal termination event"> <"Initial H3 state across locus (3/0/-1)"> <"Name tag for ouput files"> <br>
<br>
**The code generates five output text files:** <br>
(1) parameter values <br>
(2) type1a output showing fractional coverage of H3K27me3 and H3K4me1 at each timepoint, averaged over all simulated trajectories <br>
(3) type1b output showing average frequency of proximal and distal termination events at hourly timepoints, averaged over all simulated trajectories <br>
(4) type1c output showing switch-OFF times for all simulated trajectories ('0' indicates no switch-OFF) <br>
(5) type1d output showing final (last time point) spatial profile of H3K27me3 and H3K4me1, averaged over all simulated trajectories. ('1' for presence, '0' for absence) at each H3 histone <br>
<br>
**The layout of the output files are as follows:**<br>
**type1a**<br>
Column 1: Time  <br>
Column 2: H3K27me3 fractional coverage in nucleation region (averaged over all simulated trajectories) <br>
Column 3: H3K27me3 fractional coverage over the whole locus (averaged over all simulated trajectories) <br>
Column 4: H3K4me1 fractional coverage over the whole locus (averaged over all simulated trajectories) <br>
<br>
**type1b** <br>
Column 1: Time (intervals of 1h) <br>
Column 2: Number of Proximal termination events in the previous 1h (averaged over all simulated trajectories) <br>
Column 3: Number of Distal termination events in the previous 1h (averaged over all simulated trajectories) <br>
<br>
**type1c** <br>
Column 1: Switch-OFF times for all simulated trajectories <br>
<br>
**type1d** <br>
Column 1: Indices of H3 histones across the locus <br>
Column 2: H3K27me3 status at the corresponding histones (averaged over all simulated trajectories), at the last timepoint of the simulation <br>
Column 3: H3K4me1 status at the corresponding histones (averaged over all simulated trajectories), at the last timepoint of the simulation <br>
<br>
**Analysis code** <br>
The file **analysis_280423.py** contains the Python code used to analyse the simulation output <br>
<br>
**Model comparison to data** <br>
The comparison plots, such as those in Fig. 4 and Fig. 5 were made using Graphpad Prism version 10.1.0 for macOS <br>
See main text and supplementary information for details of these comparisons. <br>
Experimental datasets used in these comparisons are available on request <br> 
<br>
