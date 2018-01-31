### THIS IS AN EDITED VERSION OF THE CODE DESCRIBED BELOW. THIS IS NOT ENDORSED BY THE AUTHORS.###

====================================================
   ___ _   _ ___   _     _       ___ ___ ___ ___ 
  / __| | | |   \ /_\   | |  ___| _ ) __/ __/ __|
 | (__| |_| | |) / _ \  | |_|___| _ \ _| (_ \__ \
  \___|\___/|___/_/ \_\ |____|  |___/_| \___|___/
                                                2012
     by Jens Wetzl           (jens.wetzl@fau.de)
    and Oliver Taubmann (oliver.taubmann@fau.de)

This work is licensed under a Creative Commons
Attribution 3.0 Unported License. (CC-BY)
http://creativecommons.org/licenses/by/3.0/
====================================================

The CUDA L-BFGS library offers GPU based nonlinear
minimization implementing the L-BFGS method in CUDA.

==========================================
  USAGE
==========================================

The basic approach can be described as follows:

1. Implement your cost function in a class that
   inherits from the appropiate base class
   declared in cost_function.h

2. Create an object of class lbfgs (lbfgs.h)
   passing an object of your cost function class
   in the constructor. Adjust settings of lbfgs
   to your liking.
   
3. Run minimization providing an initial guess
   for the solution. Check the return code
   to know which stopping criterion was fulfilled.
