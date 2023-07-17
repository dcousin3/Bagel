//==============================================================================//==========//
//     _  _  _           _           _  _       _  _  _  _     _                //          //
//    (_)(_)(_)_       _(_)_       _(_)(_)_    (_)(_)(_)(_)   (_)               // comments //
//    (_)     (_)    _(_) (_)_    (_)    (_)   (_)            (_)               // the way  //
//    (_) _  _(_)   (_)     (_)   (_)  _  _    (_) _  _       (_)               // it was   //
//    (_)(_)(_)_    (_) _ _ (_)   (_) (_)(_)   (_)(_)(_)      (_)               // done in  //
//    (_)     (_)   (_)(_(_)(_)   (_)    (_)   (_)            (_)               // the old  //
//    (_)_  _ (_)   (_)     (_)   (_)_  _(_)   (_) _  _  _    (_) _  _  _       // days!    //
//    (_)(_)(_)     (_)     (_)     (_)(_)     (_)(_)(_)(_)   (_)(_)(_)(_)      //      -dc //
//                                                                              //          //
//==============================================================================//==========//
// (C) Denis Cousineau & Charles Collin, 2023                                               //
// Distributed under CC 4.0 BY NC which means:                                              //
//                BY: cite this if you use Bagel according to academic standards            //
//                NC: absolutely no commercial use without permission.                      //
//             under MIT for invFFT_1D_radix2                                               //
// Additionally, HU (humans only) which means:                                              //
//                non-humans (bots and automated web scraping algorithms) are not allowed   //
//                to use, adapt, or modify code here unless explicitely approved by the     //
//                right holders.                                                            //
//==========================================================================================//

//==========================================================================================
// Preliminary work by Anthony Liu:
//     https://github.com/turbomaze/JS-Fourier-Image-Analysis
// Massive changes by Denis Cousineau.
//     In particular, the shift was not working; and the FFT/invFFT was for 1D vectors only.
//==========================================================================================

//==========================================================================================
// Version history
//      1.2.0 (2023.05.18): Original version
//                  Detected that JS-Fourier-Image-Analysis was incorrect in many ways.
//                  For simplicity, the NxN image matrix is stored in a vector of length N*N.
//      1.3.0 (2023.05.22): Added contrast manipulation in the GUI
//                  Use (r)escaling and (m)edian equalization;
//      1.3.1 (2023.06.16): Triple-checked FFT with Mathematica 
//                  The FFT functions assumes the signal processing parameterization;
//                  Also, the matrices are coded 0..255 which has implications
//                  (relative to a 0..1 coding) when non-linear transforms are used.
//      1.3.2 (2023.06.17): Canvas are rounding pixels so we bypass them for storing images
//      1.3.3 (2023.07.02): Added truncate, rescale and removeDC
//      1.3.4 (2023.07.07): Added parameters to make letterImg
//==========================================================================================

//==========================================================================================
// Bagel is a JS library with relevant functions to perform manipulations on images.
// In particular, we have:
//      for image transforms: (all transforms returns arrays of Complex) 
//          FFT2D               a 2D fast fourier transform
//          invFFT2D            a 2D inverse FFT
//      for image generation: (all Img are arrays of Complex with real between 0 and 255)
//          emptyImg            an empty image
//          gaussianNoisyImg,   a random image
//          plaidImg            a plaid with frequencies along angles
//          gaborImg            a multiple-component Gabor patch
//          letterImg           an image of a letter
//      for image filtering: (all filters are arrays of Complex with real between 0 and 1)
//          donutFilter:        a standard band filter (low pass and high pass), 
//          gaussFilter,        a gaussian filter in the center or off the center
//          oneOverFnFilter     a 1/f^n white noise
//          bagelFilter         a bagel filter (see below) (also provided a bagel sampler)
//
// Bagel (Bubbles Along Gaussians Evolving radialLy) is a set of tools to create a filter 
//   which is composed of gaussian rings at various distances from the center. The locations 
//   are randomly chosen among the lg(width) positions with steps of 1/k with k the smooth 
//   parameter.
//   See Willenbockel et al. (2010) Journal of Experimental Psychology: Human perception & 
//        performance. doi: 10.1037/a0016465
//==========================================================================================
// Please cite: <<contact denis.cousineau@uottawa.ca for the latest reference>>.
//==========================================================================================


var Bagel = (function() {

    function version() {return "Bagel version 1.3.4"};


    //***************************************************************************************
    // image transform functions : all are programmed using the inverse FFT
    //         the results all have been tested with Mathematica
    //***************************************************************************************
    function truncate(input, lower, upper) {
        if (lower === undefined) {lower = 0;}
        if (upper === undefined) {upper = 255;}
        var out = input.map( x => ((x<lower)? lower :((x>upper)? upper :x)) );
        return out
    }
    function rescale( input, upper ) {
        // This function changes the pixel values so that they range from 0 to upper (default 255).
        // The input is assumed a square image stored in an array.
        min = Math.min.apply(null, input);
        max = Math.max.apply(null, input);
        if (min==max) {return input;}; //image is a constant so we do not touch
        if (upper === undefined) {upper = 255;}
        scaledinput = Bagel.elementWisePlus( input, -min );
        scaledinput = Bagel.elementWiseTimes( scaledinput, upper/ (max-min) );

        return scaledinput;
    }
    function recenter( input, center ) {
        // This function changes the pixel values so that its mean is the center (default 255/2).
        // The input is assumed a square image stored in an array.
        mn = Bagel.mean( input );
        if (center === undefined) {center = 255/2;}
        centeredinput = Bagel.elementWisePlus( input, center-min );

        return centeredinput;
    }
    function removeDC( input ) {
        // Ths function filters out the central, '0-frequency', aka DC, component.
        // This is only to allow better visualization of the FFT: the DC is typically
        // so large that in comparison, all the other components of the FFT are weak
        // and would show up as almost-nearly black dots. 
        // The radius is heuristically set to 1/64th the size of the whole image.
        // the input is assumed a square image stored in an array.
        dim = Math.sqrt( input.length );
        modifiedinput = Bagel.elementWiseTimes( input, Bagel.gaussFilter( [dim,dim], -dim/64 ) );
        return modifiedinput;
    }

    function FFT1D( input ) {
        // Proxy for a 1D (vector) fast-Fourier transform
        out = [];
        var tmp = input.slice();
        tmp = input.map( x => x.conjugate() );
        invFFT_1D_radix2( out, 0, tmp, tmp.length, 0, 1 )
        out = out.map( x => x.conjugate() );
        return out;
    }
    function invFFT1D( input ) {
        // Proxy for a 1D (vector) inverse fast-Fourier transform
        out = [];
        var tmp = input.slice();
        invFFT_1D_radix2( out, 0, tmp, tmp.length, 0, 1 )
        out = elementWiseTimes( out, 1/tmp.length  );
        return out;    
    }


    function FFT2D( input ) {
        // Proxy for a 2D (image) fast-Fourier transform
        // It uses the property that invFFT(x*)* is the same as FFT(x)
        // where * denotes the conjugate.
        if (input.length === 0) { throw new Error("Bagel.FFT2D:: Cannot transform an image with size of zero."); };
        if (Math.sqrt(input.length) != Math.round(Math.sqrt(input.length))) { throw new Error("Bagel.FFT2D:: Only square images supported. Current size: " + input.length); };
        if (Math.sqrt(input.length) & (Math.sqrt(input.length) - 1)) { throw new Error("Bagel.FFT2D:: Only images of dimensions power of 2 supported. Current size: " + Math.sqrt(input.length) + " x " + Math.sqrt(input.length) ); };

        out = [];
        var temp = input.slice();
        temp = temp.map( x => x.conjugate() );
        out = invFFT2D(temp)
        out = out.map( x => x.conjugate() );
        out = Bagel.elementWiseTimes( out, temp.length );

        return out
    };

    function invFFT2D( input ) {
        // Proxy for a 2D (image) inverse fast-Fourier transform
        if (input.length === 0) { throw new Error("Bagel.invFFT2D:: Cannot transform an image with size of zero."); };
        if (Math.sqrt(input.length) != Math.round(Math.sqrt(input.length))) { throw new Error("Bagel.invFFT2D:: Only square images supported. Current size: " + input.length); };
        if (Math.sqrt(input.length) & (Math.sqrt(input.length) - 1)) { throw new Error("Bagel.invFFT2D:: Only images of dimensions power of 2 supported. Current size: " + Math.sqrt(input.length) + " x " + Math.sqrt(input.length) ); };

        out = [];
        var temp = input.slice();
        invFFT_2D_radix2(out, 0, temp, [Math.sqrt(temp.length), Math.sqrt(temp.length)], 0, 1 );
        out = Bagel.elementWiseTimes( out, 1 / temp.length );
        transpose( out );

        return out
    };


    function invFFT_2D_radix2(out, start, input, dims, offset, s) {
        // It uses the property that a 2D FFT is the FFT of the column + the FFT of the rows.
        // Run invFFT_1D on rows, then on columns (by transposing the matrix)
        // With help from https://stackoverflow.com/a/76277362/5181513
        for (let i = 0; i < 2; i++) {
            for (let j = 0; j < dims[1]; j++) {
                invFFT_1D_radix2(out, j * dims[0], input, dims[0], j * dims[0], s);
            }
            input = out.slice();
            transpose( input );
        };
        out = input.slice();
    };

    function invFFT_1D_radix2(out, start, input, N, offset, s ) {
        // (c) Anthony Liu, 2014, MIT license
        if (N === 1) {
            out[start] = input[offset];
        } else {
            invFFT_1D_radix2(out, start,     input, N/2, offset,   2*s);
            invFFT_1D_radix2(out, start+N/2, input, N/2, offset+s, 2*s);
            for (let k = 0; k < N/2; k++) {
                let twiddle = cisExp(+2*Math.PI*k/N); //+
                let t = out[start+k];
                out[start+k]     = t.plus ( twiddle.times(out[start+k+N/2]) );
                out[start+k+N/2] = t.minus( twiddle.times(out[start+k+N/2]) );
            }
        }
    };


    //***************************************************************************************
    // image generation functions : emptyImg, gaussianNoisyImg, plaidImage, letterImage
    //***************************************************************************************
    function emptyImg( dims, gl ) {
        // Generate a complex array of constant pixels filled with gray level gl
        out = [];
        var N = dims[1];
        var M = dims[0];
        for (var k = 0; k < N; k++) {
            for (var l = 0; l < M; l++) {
                var idx = k*M + l;
                out[idx] = new Bagel.Complex(gl, 0);                    
            }
        };
        return out
    }

    function gaussianNoisyImg( dims, std ) {
        // Generate a complex array of zero-centered pixels coding random noise
        if (std <= 0) {throw new Error("Bagel.gaussianNoisyImg:: std must be larger than zero");};
        function gaussianRandom(mean=0, stdev=1) {
            // Standard Normal variate using Box-Muller transform.
            let u = 1 - Math.random(); // Converting [0,1) to (0,1]
            let v = Math.random();
            let z = Math.sqrt( -2.0 * Math.log( u ) ) * Math.cos( 2.0 * Math.PI * v );
            return z * stdev + mean;
        }
        out = [];
        var N = dims[1];
        var M = dims[0];
        for (var k = 0; k < N; k++) {
            for (var l = 0; l < M; l++) {
                var idx = k*M + l;
                out[idx] = new Bagel.Complex(gaussianRandom(0, std), 0);                    
            }
        };
        return out
    }

    function plaidImg( dims, angles, frequencies, radius ) {
        // Generate a complex array of gray-level pixels coding a plaid from superposition
        // of sin waves at angles with frequencies which is then given a gaussian envelop 
        // of radius.
        if (angles.length != frequencies.length) {throw new Error("Bagel.plaidImg:: There must be as many angles as frequencies.");};
        if (radius <= 0) {throw new Error("Bagel.plaidImg:: Radius must be larger than zero");};

        var plaid = Bagel.emptyImg( dims, 255/2 );
        var nplaid = 0;
        var cntr = dims[0]/2+1; //center with +1

        for (var i = 0; i < frequencies.length; i++) {
            if (frequencies[i] != 0) {
                var plaidLayer = Bagel.emptyImg( dims, 0 );
                //console.log(cntr + Math.ceil(Math.sin(angles[i]) * frequencies[i]/2-0.0001)); //caliss JS pas capable d'arrondir comme du monde...
                idx1 = (cntr + Math.ceil(Math.sin(angles[i]) * frequencies[i]/2-0.0001)) * dims[0] + ( cntr + Math.ceil(Math.cos(angles[i]) * frequencies[i]/2-0.0001) ) ;
                plaidLayer[ idx1 - 1 ] = new Bagel.Complex(1, 0);
                idx2 = (cntr - Math.floor(Math.sin(angles[i]) * frequencies[i]/2+0.0001)) * dims[0] + ( cntr - Math.floor(Math.cos(angles[i]) * frequencies[i]/2+0.0001) ) ;
                plaidLayer[ idx2 -1 ] = new Bagel.Complex(1, 0);
                //console.log(idx1,idx2);
                plaidLayer = Bagel.invFFT2D( plaidLayer );
                plaidLayer = plaidLayer.map( x => x.magnitude() ); //remove phase information
                plaidLayer = Bagel.rescale( plaidLayer );  // onto 0..255
                plaid = Bagel.elementWisePlus(plaid, plaidLayer);
                nplaid ++;
            }
        };
        if (nplaid > 0) { plaid = Bagel.elementWiseTimes( plaid, 1/(nplaid+0.5)); }; //0..255
        plaid  = Bagel.elementWisePlus( plaid, -255/2 );                             //-255/2 .. +255/2
        filter = Bagel.gaussFilter( dims, radius );
        plaid  = Bagel.elementWiseTimes( plaid, filter );
        plaid  = Bagel.elementWisePlus( plaid, +255/2 );                             //0..255

        // not symmetrically located??
        mn = Bagel.min(plaid.map(x => x.magnitude() ) );
        mx = Bagel.max(plaid.map(x => x.magnitude() ) );
        if (nplaid > 0) { var tt = (255-mx)/2 - (mn-0)/2;
        plaid  = Bagel.elementWisePlus( plaid, tt ); }

        return plaid;    
    }

    function gaborImg( dims, angles, frequencies, radius ) {
        // Generate a complex array of gray-level pixels coding a gabor from superposition
        // of sin waves at angles with frequencies which is then given a gaussian envelop 
        // of radius.
        if (angles.length != frequencies.length) {throw new Error("Bagel.plaidImg:: There must be as many angles as frequencies.");};
        if (radius <= 0) {throw new Error("Bagel.gaborImg:: Radius must be larger than zero");};

        var gabor = Bagel.emptyImg( dims, 0 );
        var cntr = dims[0]/2+1; //center  +++ VALIDATED WITH math

        for (var i = 0; i < frequencies.length; i++) {
            if (frequencies[i] != 0) {
                idx1 = (cntr + Math.round(Math.sin(angles[i]) * frequencies[i]/2)) * dims[0] + ( cntr + Math.round(Math.cos(angles[i]) * frequencies[i]/2) ) ;
                gabor[ idx1 - 1 ] = new Bagel.Complex(1, 0);
                idx2 = (cntr - Math.round(Math.sin(angles[i]) * frequencies[i]/2)) * dims[0] + ( cntr - Math.round(Math.cos(angles[i]) * frequencies[i]/2) ) ;
                gabor[ idx2 -1 ] = new Bagel.Complex(1, 0);
            }
        };
        // inverting the points
        gabor = Bagel.invFFT2D( Bagel.shiftQuadrants( gabor ) );
        temp = gabor.map( x => x.magnitude() ); //remove phase information
        tempMax = Math.max.apply(null,temp);
        gabor = gabor.map( x => new Bagel.Complex(x.real* 255/tempMax, x.imag * 255/tempMax) );

        // putting the gaussian envelop
        gabor  = Bagel.elementWisePlus( gabor, -255/2 );
        filter = Bagel.gaussFilter( dims, radius );
        gabor  = Bagel.elementWiseTimes( gabor, filter );
        gabor  = Bagel.elementWisePlus( gabor, +255/2 );

        return gabor;    
    }

    function letterImg( dims, letter, txtformat, bgcolor ) {
        // Generate a complex array of b/w pixels coding a given letter.
        // Uses an in-memory canvas for its ability to draw text
        if (txtformat === undefined) {txtformat = "bold 128px Arial"; };
        if (bgcolor === undefined) {bgcolor = 255; }; // actually only gray level from 0:black to 255:white

        // locate a font dimension in px in the format
        var size = txtformat.split(/\s/).map( x => parseInt(x, 10) ).filter( x => x)[0];
        //console.log("located font size: ", size );

        var inMemoryCanvas = document.createElement('canvas');
        inMemoryCanvas.width  = dims[0];
        inMemoryCanvas.height = dims[1];
        var imc = inMemoryCanvas.getContext('2d'); 
        imc.textBaseline = "top";
        imc.font      = txtformat;              // change as desired
        imc.clearRect(0,0, dims[0], dims[1] );  // erase all
        imc.fillText(letter, (dims[1] - imc.measureText(letter).width)/2, (dims[0]-0.75*size)/2 );
        var letterImageData = imc.getImageData(0, 0, dims[0], dims[1] ).data;
        // remove color
        var bwletterImageData = [];
        for (var i = 0; i < letterImageData.length; i+=4) {
            // it is bizzarely coded with the transparent layer only... the letters are always black
            bwletterImageData[i/4] = Math.min(bgcolor, 255-letterImageData[i+3] );
        }
        // convert to complex
        bwletterImageData = bwletterImageData.map( x => new Bagel.Complex(x, 0) )  ;
        return bwletterImageData
    }


    //***************************************************************************************
    // image filter functions : donut, gauss, oneOverFn, bagel
    //***************************************************************************************
    function donutFilter( dims, lowPass, highPass ) {
        // AKA band-pass filter: a hard-threshold filter inside lowPass-highPass
        out = [];
        // erase the pixels that are within the threshold distances
        var lowPassSq = Math.pow(lowPass, 2);
        var highPassSq = Math.pow(highPass, 2);
        var N = dims[1];
        var M = dims[0];
        for (var k = 0; k < N; k++) {
            for (var l = 0; l < M; l++) {
                var idx = k*M + l;
                var d = Math.pow(k-M/2, 2) + Math.pow(l-N/2, 2);
                out[idx] = new Bagel.Complex(0, 0);                    
                if ( // low pass only (e.g., ",4") or high pass only (e.g., "4,")
                    (d >= lowPassSq && isNaN(highPass)  ) ||
                    (d <= highPassSq && isNaN(lowPass)  ) 
                ) { out[idx] = new Bagel.Complex(1, 0); }
                if ( // band pass filter (e.g, "4,12")
                    d >= lowPassSq && d <= highPassSq  && !isNaN(lowPass) && !isNaN(highPass) 
                ) { out[idx] = new Bagel.Complex(1, 0); }
            }
        };
        return out
    }

     function gaussFilter( dims, sigma ) {
        // A soft filter: a gaussian function centered on zero with radius sigma
        // Use a negative sigma to have the complement.
        out = [];
        radius = Math.abs(sigma);
        inout  = Math.sign(sigma);  // +1 => 0*255 + 1*C; -1 => 1*255 -1*C
        for (var k = -dims[0]/2; k < dims[0]/2; k++) {
            for (var l = -dims[1]/2; l < dims[1]/2; l++) {
                reel =  round(Math.exp(-(Math.pow(k,2)+Math.pow(l,2) )/ Math.pow(radius,2)),14);
                out[ (k+dims[0]/2)* dims[0] + l+dims[0]/2 ] = 
                    new Bagel.Complex( (1-inout)/2 * 1 + inout * reel, 0 )
            }
        }
        return out
    };

     function oneOverFnFilter( dims, n ) {
        // 1/f^n filter: an exponentially-decreasing function centered on zero with power n.
        out = [];
        for (var k = -dims[0]/2; k < dims[0]/2; k++) {
            for (var l = -dims[1]/2; l < dims[1]/2; l++) {
                d = Math.sqrt( Math.pow(k,2)+Math.pow(l,2) )
                reel = round( 1 / Math.pow(d, n ), 6);
                out[ (k+dims[0]/2)* dims[0] + l+dims[0]/2 ] = 
                    new Bagel.Complex( Math.min( 1, reel), 0 )
            }
        }
        return out
    };

    function bagelSampler( dims, b, k ) {
        // Sample b frequencies from the range 0..2w, where w is half the size dims[0], 
        // in increment of k (the smoothing parameter) using a logarithmic scale based on octaves.
        if (b <= 0) {throw new Error("Bagel.bagelSampler:: Sample size b must be larger than zero");};
        if (k <= 0) {throw new Error("Bagel.bagelSampler:: Smooth parameter k must be larger than zero");};

        w = dims[0]/2;  // half width of the image
        freqs = [];     // a vector with the b frequencies sampled
        indicator = []; // a vector of 0s with 1s where a freq has been sampled
        temp = [];
        for (var i = 0; i <= k*Math.log2(2*w); i++) { 
            temp.push(round(i,3)); 
            indicator.push(0);
        };
        temp = shuffle(temp);
        for (var i = 0; i < b; i++) {
            freqs.push( 2**(temp[i]/k) ) ;
            indicator[temp[i]]=1;
        };
        freqs.sort( (a,b) => a-b );
        return [freqs, indicator];
    };

    function bagelFilter( dims, mus, sigmas ) {
        // Generate a bagel filter, that is, radial gaussians located on a certain frequency
        // (a certain distance from the center of the filter). The lists mus and sigmas must
        // have the same length. It is suggested to have the sigmas = mus / (2*k).
        // The bagels are centered on (129,129) if dims = [256,256] 
        // but 129 is the 128th column (as we start from zero) 
        if (mus.length != sigmas.length) {throw new Error("Bagel.bagelFilter:: There must be as many mus as sigmas.");};

        var w    = dims[0]/2;
        var out  = new Array(dims[0]*dims[1]);   
        for (let i = 0; i < dims[0]*dims[1]; i++) out[i] = 0;

        // compute the sum of the gaussian distributions along a line of length 2w
        var curve = []; //in the center, curve heigth is null
        for (let x = 0; x < 2*w; x++) {
            var h = 0; 
            for (let j = 0; j < mus.length; j++ ) {
                h += Math.exp(-1 * (x-mus[j])**2 / (2* sigmas[j]**2 ) );
            }
            curve.push(h)
        }
        curve = curve.map( x => round(x, 10))

        // rotate this mixture around the center of one quadrant; 
        // it is symmetrical about the other three quadrants
        for (let k = 0; k < 2*w; k++) {
            for (let l = 0; l < 2*w; l++) {
                var d = Math.round( Math.sqrt( (k-w)**2 + (l-w)**2 ) );
                out[ k*dims[1] + l ] = curve[Math.max(d,1)];
            }
        };
        out = rescale(out, 1); // top the filter to 1
        out = out.map( x => new Complex( x, 0 ) );
        return out;
    };



    //***************************************************************************************
    // shift quadrants function.
    //***************************************************************************************
    function shiftQuadrants( matrix ) {
        // Shift the quadrants of the matrix.
        // Note that shift of shift returns the orignal matrix. 
        var out = [...matrix]; 
        var dims = [Math.sqrt(matrix.length), Math.sqrt(matrix.length)];
        var M = dims[0];
        var N = dims[1];

        for (var n = 0; n < N/2; n++) {
            for (var m = 0; m < M/2; m++) {
                // swap upper left quadrant
                out[(n+0)*N   + (m+0)  ] = matrix[(n+N/2)*N + (m+M/2)]
                // swap lower right quadrant
                out[(n+N/2)*N + (m+N/2)] = matrix[(n+0)*N   + (m+0)  ]
                // swap upper right quadrand
                out[(n+0)*N   + (m+N/2)] = matrix[(n+N/2)*N + (m+0)  ]
                // swap lower left quadrant
                out[(n+N/2)*N + (m+0)  ] = matrix[(n+0)*N   + (m+M/2)]
            }
        }
        return out
    };


    //***************************************************************************************
    // Subsidiary functions:
    //     (math) cisExp, median, round, (matrix) transpose, (random) shuffle
    //***************************************************************************************
    function cisExp(x) { // e^ix = cos x + i*sin x
        return new Complex(Math.cos(x), Math.sin(x));
    };
    function round(n, places) {
        var mult = Math.pow(10, places);
        if (Math.abs(n) < 1/mult) {
            return 0.00
        } else {
            return Math.round(mult*n)/mult;
        }
    };
    function min(numbers) {
        return Math.min.apply(null, numbers)
    }
    function max(numbers) {
        return Math.max.apply(null, numbers)
    }
    function mean(numbers) {  //  inspired by https://stackoverflow.com/a/62372003/5181513
        return numbers.reduce((acc,x)=>(acc+x), 0) / numbers.length;
    }
    function median(numbers) {
        // see https://stackoverflow.com/a/53660837/5181513
        if(numbers.length === 0) throw new Error("Bagel.median:: No inputs to compute median");
        
        const sorted = Array.from(numbers).sort( (a, b) => a - b );
        const middle = Math.floor(sorted.length / 2);
        if (sorted.length % 2 === 0) {
            return (sorted[middle - 1] + sorted[middle]) / 2;
        } else {
            return sorted[middle];
        }
    };
    function transpose( matrix ) {
        // matrix is a square matrix kept in a vector
        dims = [ Math.sqrt(matrix.length), Math.sqrt(matrix.length) ]
        const swap = (arr, i, j) => {
            const t = arr[i];
            arr[i] = arr[j];
            arr[j] = t;
        };
        for (let i = 1; i < dims[1]; i++) {
            for (let j = 0; j < i; j++) {
                swap(matrix, i * dims[0] + j, j * dims[0] + i);
            }
        }
    };
    function shuffle(array) {
        // Shuffle the array randomly; from https://stackoverflow.com/questions/2450954
        var out = array.slice()
        var currentIndex = out.length, temporaryValue, randomIndex ;
        // While there remain elements to shuffle...
        while (0 !== currentIndex) {
            // Pick a remaining element...
            randomIndex = Math.floor(Math.random() * currentIndex);
            currentIndex -= 1;
            // And swap it with the current element.
            temporaryValue      = out[currentIndex];
            out[currentIndex] = out[randomIndex];
            out[randomIndex]  = temporaryValue;
        }
        return out;
    }
    function imageDiags( image, label ) {
        if (image.constructor.name === "Array") {
            var dims = Math.sqrt( image.length );
            var temp = image;
            if (image[0].constructor.name === "Complex") { 
                temp = temp.map( x=>x.magnitude() );
            };
            console.log( "   ", label, ": ");
            console.log( "      => type      = ", image[0].constructor.name);
            console.log( "      => Dimension = ", dims);
            console.log( "      => Min, med/mean max  = ", [Math.min.apply(null, temp), Bagel.median(temp), Bagel.mean(temp), Math.max.apply(null, temp)] );
        }
    }


    //***************************************************************************************
    // Basic vector operators: .+, .*, .=
    //     These functions matches entries of an array to entries of a second array if
    //     it is an array of the same type, or to the second argument if it is a single value.
    //***************************************************************************************
    function elementWisePlus(a,b){
        if ((a.constructor.name != "Array")) {
            throw new Error("Bagel.elementWisePlus: First argument must be an array...");
        };
        if (a[0].constructor.name === "Complex") {
            if (b.constructor.name === "Number") {
                return a.map( (x,i) => x.plus( b ) );
            } else if (b.constructor.name === "Array") {
                return a.map( (x,i) => x.plus( b[i] ) );
            } else {
                throw new Error("Bagel.elementWisePlus: Second argument must be a scalar or an array...");
            }
        } else if (a[0].constructor.name === "Number") {
            if (b.constructor.name === "Number") {
                return a.map( (x,i) => x + b  );
            } else if (b.constructor.name === "Array") {
                return a.map( (x,i) => x + b[i]  );
            } else {
                throw new Error("Bagel.elementWisePlus: Second argument must be a scalar or an array...");
            }
        }
    };
    function elementWiseTimes(a,b){
        if ((a.constructor.name != "Array")) {
            throw new Error("Bagel.elementWiseTimes: First argument must be an array...");
        };
        if (a[0].constructor.name === "Complex") {
            if (b.constructor.name === "Number") {
                return a.map( (x,i) => x.times( b ) );
            } else if (b.constructor.name === "Array") {
                return a.map((x,i) => x.times( b[i] ) );
            } else {
                throw new Error("Bagel.elementWiseTimes: Second argument must be a scalar or an array...");
            } 
        } else if (a[0].constructor.name === "Number") {
            if (b.constructor.name === "Number") {
                return a.map( (x,i) => x * b  );
            } else if (b.constructor.name === "Array") {
                return a.map((x,i) => x * b[i] );
            } else {
                throw new Error("Bagel.elementWiseTimes: Second argument must be a scalar or an array...");
            } 
        }
    };
    function elementWiseEqual(a,b){
        if ((a.constructor.name != "Array")) {
            throw new Error("Bagel.elementWiseEqual: First argument must be an array...");
        };
        if (a[0].constructor.name === "Number") {
            if (b.constructor.name === "Number") {
                return a.every( (x,i) => x === b );
            } else if ((b.constructor.name === "Array")&&(b[0].constructor.name === "Number")) {
                return a.every( (x,i) => x === b[i])
            } else {
                throw new Error("Bagel.elementWiseEqual: Second argument must be a scalar or an array of scalar to match the first argument...");
            }
        } else if (a[0].constructor.name === "Complex") {
            if (b.constructor.name === "Complex") {
                return (a.every( (x,i) => x.real === b.real ) && a.every( (x,i) => x.imag === b.imag ));
            } else if ((b.constructor.name === "Array")&&(b[0].constructor.name === "Complex")) {
                return (a.every( (x,i) => x.real === b[i].real) && a.every( (x,i) => x.imag === b[i].imag));
            } else {
                throw new Error("Bagel.elementWiseEqual: Second argument must be a complex or an array of complex to match the first argument...");
            }
        }
    };
    

    //***************************************************************************************
    // data structure Complex.
    //     The Complex class implements complex numbers along with some basic operations.
    //***************************************************************************************
    class Complex{
        constructor(re, im) {
            this.real = re;
            this.imag = im;
        };
        magnitude = function() {
            return Math.sqrt(this.real*this.real + this.imag*this.imag);
        };
        conjugate = function() {
            return new Complex( this.real, -this.imag );
        };
        plus = function(z) {
            if ( z.constructor.name === "Complex" ) { // complex multiplication
                return new Complex( this.real+z.real, this.imag+z.imag );
            } else { // scalar plus
                return new Complex( this.real+z, this.imag);
            }
        };
        minus = function(z) {
            if ( z.constructor.name === "Complex" ) { // complex multiplication
                return new Complex( this.real-z.real, this.imag-z.imag );
            } else { // scalar minus
                return new Complex( this.real-z, this.imag );
            }
        };
        times = function(z) {
            if ( z.constructor.name === "Complex" ) { // complex multiplication
                var rePart = this.real*z.real - this.imag*z.imag;
                var imPart = this.real*z.imag + this.imag*z.real;
                return new Complex(rePart, imPart);
            } else { // scalar multiplication
                return new Complex( z*this.real, z*this.imag);
            }
        };
    }


    //***************************************************************************************
    // returned value 
    //***************************************************************************************
    return {
        version:            version,        // returns the version of Bagel
        Complex:            Complex,        // class to store complex number

        round:              round,          // round number to a certain number of decimals
        min:                min,            // shorter than Math.min.apply(null, image)
        max:                max,
        median:             median,         // compute the median of an array
        mean:               mean,
        transpose:          transpose,      // transpose a square maxtrix stored in an array
        shuffle:            shuffle,        // randomize the items' order in an array
        imageDiags:         imageDiags,     // some basic diagnostic on the image content

        truncate:           truncate,       // truncate entries of a real matrix between lower and upper (default 0 and 255)
        rescale:            rescale,        // rescale between 0 and 255 the items from an array
        recenter:           recenter,       // move the pixels so that the mean is center (default 255/2)
        removeDC:           removeDC,       // fitler out the center of a FFT image
        FFT1D:              FFT1D,          // perform fast-Fourier transfomr on a vector
        invFFT1D:           invFFT1D,       // perform an inverse fast-Fourier transform on a vector
        FFT2D:              FFT2D,          // perform fast-Fourier transform on a matrix
        invFFT2D:           invFFT2D,       // perform an inverse fast-Fourier transform on a matrix
        shiftQuadrants:     shiftQuadrants, // shift the quadrants from a matrix

        donutFilter:        donutFilter,    // filter composed of a band of ON (1) on a matrix of OFF (0)
        gaussFilter:        gaussFilter,    // filter composed of a mode culminating at 1 in the center of the image
        oneOverFnFilter:    oneOverFnFilter, // fiter compose of an exponential with exponent -n over the center of the image
        bagelFilter:        bagelFilter,    // filter composed of a ring whose transversal shape is a gaussian
        bagelSampler:       bagelSampler,   // adjuvent to bagelFilter which picks bagels at random

        gaussianNoisyImg:   gaussianNoisyImg, // Normally-distributed random gray-level pixels
        emptyImg:           emptyImg,       // empty image...
        letterImg:          letterImg,      // pixels representing a letter in black on a white background
        plaidImg:           plaidImg,       // an image obtained from a superposition of one-frequency images
        gaborImg:           gaborImg,       // an image obtained from the inverse FFT of multiple frequency components

        elementWisePlus:    elementWisePlus,    // add items from two matrix or one matrix and one scalar
        elementWiseTimes:   elementWiseTimes,   // mutiply items ...
        elementWiseEqual:   elementWiseEqual    // check equality of items ...
        
    };
})();

console.log( Bagel.version() );
