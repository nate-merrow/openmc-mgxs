## MGXS Project
I created this project to gain an understanding of how to generate multigroup cross sections in openmc.

The main script (written in a python file through jupyter notebook), openmc-mgxs.py, imports the openmc program, 
defines materials, geometry, settings, tallies, and the multigroup cross section library, and then generates several plots using matplotlib.
One is a heat map of the scattering matrix and the other three are bar charts of the total, fission, and absorption cross sections [in barns] of the fuel.
There is also the tallies.out file which contains some numerical output defined in the tallies section of the python script

For this project, I defined three neutron energy groups based on E.E.Lewis's "Fundamentals of Nuclear Reactor Physics", including  fast, intermediate, and thermal spectra.
Ideally, the cross-sectional data from this openmc project would then be exported to an hdf5 file that would then be used as an input file for a multiphysics simulation engine like MOOSE, 
in which you could then do a variety of diffusion, transport, performance, aging, etc. simulations.

## From Here
I hope to install MOOSE in the future to potentially see this type of project all the way through.
I may also implement similar mgxs code in a mox fuel project with similar energy groups to see how the fuels cross sections vary by energy group.
