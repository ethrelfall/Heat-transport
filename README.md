# Heat-transport
Models of heat transport

**Ed_heat_transport_toy.py** - toy model of heat transport given a temperature profile on the left-hand-side of the domain, insulating boundary conditions on the left-hand-side half of the top and bottom, then constant-temperature boundary conditions for the right-hand-side portion of the domain.

The diffusivity tensor is anisotropic on the LHS (with transport only permitted in the horizontal direction, or at 5 degrees to that axis, in the examples illustrated).  Observe, however, how contact with the RHS region "sucks" heat out of the beam, reducing the temperature at the strike point compared to the LHS wall temperature.

![heat_transport_toy_output](png/Ed_heat_transport_toy_output.png "Output of heat transport toy for anisotropic (50,0) on LHS and isotropic (1,1) on RHS.")

![heat_transport_toy_output_5deg](png/Ed_heat_transport_toy_output_5deg.png "Output of heat transport toy as above, with anisotropic diffusion axis aligned at 5 deg to the horizontal (note the LHS temperature profile has been adjusted to a baseline of zero to avoid boundary artifacts)).")

It is easy to create more complicated examples e.g. anisotropy tensor aligned to curved field lines (e.g. change the definition of bhat to bhat = as_vector([1.0, sin(20.0*x)])) to create the following plot.

![heat_transport_toy_output_curved](png/Ed_heat_transport_toy_output__curved.png "Output of heat transport toy with anisotropy aligned to curved field lines.")

Some of the basic physics can be shown in a 1D example (Neumann BCs top and bottom) with unit temperature on the left-hand-side boundary and zero on the right.  Then the temperature varies piecewise linearly and the midpoint temperature is k_left / (k_left+k_right) i.e. the fraction of temperature dropped across the plasma depends on the ratio of its conductivity to that of the metal, and the more conductive the plasma the less temperature is dropped across it. Here k_left is 5 times k_right and so the strike-point temperature is 5/6 = 0.833.  Note the heat flux is constant in x as it is the product of k and the gradient (factors which vary reciprocally).

![heat_transport_toy_output_1d](png/Ed_heat_transport_toy_output_1d.png "Output of heat transport toy for 1D scenario with k=5 on the left half and 1 on the right.")
