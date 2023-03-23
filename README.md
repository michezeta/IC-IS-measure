# IC-IS-measure
Codes to reproduce the IC-IS measure I proposed in Zema (2023) for price discovery in a non-Normal setting.

The folder contains the following files:

1) dataimport&processing: Used to process the high-frequency data to get a final clean dataset on which the IC-IS should be estimated for the empirical application. You can use this script to process your event-time or natural-time based prices data. Notice that the data behind the empirical application cannot be shared.

2) set.of.procedures: It contains a set of fundamental functions needed to perform the PML estimates of the mixing matric C.

3) main: Main script used to perform the empirical analysis and estimates the IC-IS on empirical data. This script calls 'set.of.procedures' to perform the estimates.

4) simulations: Script containing the code used to perform the simulation section of the paper. Also this script calls 'set.of.procedures' to perform the estimates on the simulated data.

5) useful functions: Additional script containing some useful (but not necessary) functions. 
