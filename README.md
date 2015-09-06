Small collection of codes for running a simulation with a discrete-time sliding mode controller.
This repository features two implementations: one using MATLAB and the other one in Python.

# MATLAB code

The MATLAB files can be used to implement such a controller in Simulink inside an ``embedded MATLAB function''.
It is also possible to just run a simulation with them.

For now we only support the twisting algorithm and a classical SMC controller with diagonal CB matrix.

# Python code

The Python code extensively uses the bindings of the [Siconos](http://github.com/siconos/siconos) platform.
You need to install this dependency before being able to run the code. The files presented here make use of
the Control module of Siconos. It would be possible to directly formulate the control input value as the
solution of an Affine Variational Inequality (AVI) and use the solver from the Numerics module.

The supported controller are:
- classical SMC (with or without equivalent part)
- twisting controller (with or without the modification to make it finite-time stable)
