# Heat-transport
Models of heat transport

**Ed_heat_transport_toy.py** - toy model of heat transport given a temperature profile on the left-hand-side of the domain, insulating boundary conditions on the left-hand-side half of the top and bottom, then constant-temperature boundary conditions for the right-hand-side portion of the domain.

The diffusivity tensor is anisotropic on the LHS (with transport only permitted in the horizontal direction, or at 5 degrees to that axis, in the examples illustrated).  Observe, however, how contact with the RHS region "sucks" heat out of the beam, reducing the temperature at the strike point compared to the LHS wall temperature.

![heat_transport_toy_outout](png/Ed_heat_transport_toy_output.png "Output of heat transport toy for anisotropic (50,0) on LHS and isotropic (1,1) on RHS.")

![heat_transport_toy_outout_5deg](png/Ed_heat_transport_toy_output_5deg.png "Output of heat transport toy as above, with anisotropic diffusion axis aligned at 5 deg to the horizontal (note the LHS temperature profile has been adjusted to a baseline of zero to avoid boundary artifacts)).")
