# Sphere-Steady-Flow
MATLAB code that solves and compiles the steady flow around a rotating sphere. The package contains code that solves the rotating sphere boundary layer, a new model of the impinging region and the subsequent eruption of the radial jet.

The impinging region has been rederived via a new, but rudimentary, scaling argument but is still solved in the same manner as the [Stewartson](https://github.com/bensmith95/Stewartson-1958-Sphere-Flow) code, via a multigrid method. However, the azimuthal momentum equation is now coupled which does slow convergance slightly, nevertheless, it still is significantly quicker than a full DNS. Additionally, the radial jet is also solved via both Picard and Newton methods; the Picard method is needed to obtain a reasonable guess for the Newton method to find a convergent solution. For full details see Smith et al (2022) [in draft].

To run the code, download all the files making sure the directory structure remains intact. Then simply enter _BaseFlow(**Re**)_ in the command window, where _**Re**_ is the square root of the Reynolds number. The data will be saved in a directory called _Flows_ that will be created automatically. To see the results run the file _figures.m_.

<p align="center">
  <img width="650" src="https://user-images.githubusercontent.com/29705711/177974288-59705a65-884c-471b-b820-ce32bc1dccd2.png">
</p>
