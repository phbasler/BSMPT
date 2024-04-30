In case you want to use certain classes/functions of BSMPT without using the whole program we provide some examples on how to do it.
To use it just put your `.cpp` file inside this folder and (re-)compile BSMPT.

The executables will be placed in `bin/standalone/` inside the build-folder.

We provide a few examples:
- `CalculateAction.cpp` - Solves the bounce equation and calculate the Euclidian action. The user is expected to provide the initial guess path and the potential, the gradient is optional.
- `GenericModel.cpp` - The user provides a potential \f$V(\phi)\f$, the zero-temperature VEV and the dimensionality of the VEV. This tracks the minima and calculates characteristic temperatures as well as the GW spectrum of first-order phase transitions.
- `TunnelingPath.cpp` - Solve the bounce equation using the full `BSMPTv3` and prints the tunneling path and the VEV profile in `Mathematica` and `python` formats.
