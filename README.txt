############################################################################################
				NINE Propagator
				
				version 5 Kinetic Impactor / specially modified for Itokawa spin up 
				
				by Siegfried Eggl 2012-2015
					
############################################################################################


!!!!WARNING!!!!: 	

This is not an "official release"! 
The code is under development and may contain any number of nasty bugs, especially if the input is messed with!
Also, it is in no way "clean" i.e. not all is commented and some programming solutions are fairly horrible...
Should you encounter any difficulties I shall be happy to assisst, however.

One more thing:The intention of the current version is post-mitigation risk analysis. 
Therefore, most of the code contains Kinetic Impactor passages that should not be tempered with. 
This function can be switched off in the config.inn file, though, by setting the "kick" variable to "n".

Be aware that the close encounter initial date is still active!
					
					
					
COMPILATION:

The source file modules are in the "src_new" folder, a make file is in the current directory.
currently gfortran and ifort compilers are supported with -O optimization. 
For some reason -O2 and -O3 optimizations do no longer work on the latest versions of gcc and ifort?!?!?

-> go into the nine_prop_v4 directory: type "make clean" -> gets rid of all precomiled modules;
				       type "make"       -> the program "nine.exe" should magically appear


INPUT:

FOUR things are required: 

1) JPL ephemeris (binary, linux/unix, version >= DE405). The file has to be called "JPLEPH" and it has to rest in the "ephem" folder.
2) "config.inn" file: this file contains all sorts of parameters for the integration, i.e. start and stop dates, precision, etc. A sample file is provided,
it has to be in the same directory as the "nine.exe" executable
3) "massive.inn" file: this file contains the initial conditions of massive asteroid perturbers that interact gravitationally among each other and with the massless particles.
Also this file has to be in the same directory as "nine.exe"
4) Finally, the "massless.inn" file is basically the same as the massive.inn, only for massless particle, i.e. no masses are required.

Generally, all units are Gaussian, i.e. au, D, k, ...

If different initial epochs are given for the massive and massless particles, the integrator propagates everything to "tstart" given in the config.inn file.
However, backwards integration does not work yet! Sorry. So you better make sure that tstart is > than the initial epochs given in the massive and massless files.

EXECUTION
 
If all input files are present: "./nine.exe" should start the calculations on Unix machines.


OUTPUT:

Some different output files can be chosen in the "config.inn" file. A list of the available files is given in the config.inn file.
They are selected by adding keywords to the last line in the config.inn file (see sample config.inn for details)

Generally there can be output for ephemeris, massive and massless particles.

The files that are most interesting for impact risk analysis are the ".dmi" files that are going to be produced if "cen" is selected as an output, where the minimum encounter distances with Earth are recorded (5th and 8th column).
The 5th column is the more reliable interpolation via splines, the 8th column is usually too pessimistic (i.e. too small encounter distances) since it only uses
initial conditions at the sphere of influence and not boundary conditions as splines do.
Encounter velocities and dates are given as well.

The b-plane coordinates do not work - sorry, did not get around to fix them yet.

IMPACT PROBABILITIES

If you want to calculate impact probabilities, you need to associate a priori uncertainties to your colones in the massless.inn file.
I did that using ORBFIT to generate clones along the LOV. They provide the respective sigma (gauss distribuition) to every single clone in the ".fou" file, so that
you can estimate the impact probability by 

P=Sum_i [1/2 * { - ERF(sigma_i/sqrt(2)) + ERF(sigma_i+1/sqrt(2))}], 

where ERF is the Error function, and i and i+1 are two consecutive clones that are below the impact minimum distance (e.g. 1 R_earth)

Sample programs of how to create a "massless.inn" file from an Orbfit ".fou" output and how to extract impact probabilites can be found in the "additional" folder.  




BUGS FIXED:
18.05.2016 JD to Gregorian date conversion
08.02.2015 now works with DE431 ephemeris and logarithmic delta-V sampling if the difference between imparted velocities is more than one order of magnitude
06.08.2013 close encounter and deep close encounter numbers were wrong 



Have fun!





