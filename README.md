Bagel: JavaScript library for image manipulation
==

This library has two purposes, (a) to provide functions to perform image manipulation;
(b) to run Bagel image generation.

Bagel (Bubbles Along Gaussians Evolving radialLy) is a set of tools to create a filter 
which is composed of gaussian rings at various distances from the center. The locations 
are randomly chosen among the lg(width) positions with steps of 1/k with k the smooth 
parameter.

Preliminary work were presented in Willenbock et al. (2010). This set of tools expands on
this work by providing functional JS code, and by offering a slightly different sampler
based on a log_2 space of frequencies.
   
You can find a live demo page at: http://dcousin3.github.io/Bagel/Bagel.htm This page lets you
test various functions using a rudimentary interface.

   
## License
Bagel is distributed under a CC 4.0 BY NC license
One exception is the function `invFFT_1D_radix2` which is distributed under a MIT license.

 
## Reference

Cousineau, D., & Collin, C. (2023) Bagel: JavaScript library for image manipulation [software]
    version 1.3.1. URL: http://dcousin3.github.io/Bagel/

Willenbockel et al. (2010) Journal of Experimental Psychology: Human perception & 
		performance. doi: 10.1037/a0016465 