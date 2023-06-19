//************************************************************
//
// Implementing functions that returns bagels
// (2023) Cousineau & Collins
// Work under way. Comments welcome.    
//
//************************************************************


var BagelImageAnalysis = (function() {
    //-----------------------------------------------------
    // global variables 
    //-----------------------------------------------------
    var canvases = [];  // canvases on the html page
    var ctxs     = [];  // contexts within the canvases

    var dims     = [-1, -1]; // will be set later

    var rawImage = [];  // (1) original image
    var newImage = [];  // (2) transformed image 
    var cstImage = [];  // (3) contrasted image
    var filter   = [];  // a filter

    
    //-----------------------------------------------------
    // initializing the listeners for all the buttons
    //-----------------------------------------------------
    function initBagelImageAnalysis() {

        $s('#draw-cs-btn').addEventListener('click', function() {
            disableButtons(loadImage('https://raw.githubusercontent.com/turbomaze/JS-Fourier-Image-Analysis/master/cs.png'));
        });

        $s('#draw-grace-btn').addEventListener('click', function() {
            disableButtons(loadImage('https://raw.githubusercontent.com/turbomaze/JS-Fourier-Image-Analysis/master/grace.png'));
        });

        $s('#draw-lettre-btn').addEventListener('click', function(){
            disableButtons(generateLetter(lettre.value));
        });

        $s('#draw-plaid-btn').addEventListener('click', function(){
            disableButtons(generatePlaid(plaidfrq.value));
        });

        $s('#transform-btn').addEventListener('click', function() {
            if (dims[0]==-1) { return alert('You need to draw an image to canvas 1 first.'); }
            disableButtons(transformAction(arg.value, arglist.value, noiseprop.value, nselist.value ));
        });

        $s('#method-btn').addEventListener('click', function() {
            if (dims[0]==-1) { return alert('You need to transform an image to canvas 2 first.'); }
            disableButtons(contrastAction(method.value, range.value ));
        });

        // initialize the canvas found on the html document
        for (var ai = 0; ai < 3; ai++) {
            canvases[ai] = $s('#canvas'+ai);
            ctxs[ai] = canvases[ai].getContext('2d');
        }
        console.log("End of initializations. Listening...");

    }


    //-----------------------------------------------------
    // two shortcut functions
    //-----------------------------------------------------
    function adjustCanvases( dims ) {
        for (var ai = 0; ai < 3; ai++) {
            canvases[ai].width = dims[0];
            canvases[ai].height = dims[1];
        }
    }
    function someDiags( image ) {
        console.log( "      => Dimension = ", dims);
        console.log( "      => Min, max  = ", [Math.min.apply(null, image), Math.max.apply(null, image)] );
    }
    function fillCanvasFrom( ctx, matrix ) {
        ctx.clearRect(0, 0, 256, 256); // erase all
        var ImageData = ctx.getImageData(0, 0, dims[0], dims[1] );
        for (var i = 0; i < dims[0]; i++) {
            for (var j = 0; j < dims[1]; j++) {
                var idx = i*dims[1] + j;
                for (var c = 0; c < 3; c++) {
                    ImageData.data[4*idx+c] = matrix[idx] ;
                }
                // full alpha on the fourth channel
                ImageData.data[4*idx+3] = 255 ;
            }
        }
        ctx.putImageData(ImageData, 0, 0);
    }


    //-----------------------------------------------------
    // load or generate to rawImage and show in canvas0
    //-----------------------------------------------------
    function loadImage( loc ) {
        console.log("===Step 1a: Loading an image===");

        // create an in-memory image which will fetch the URL's image
        var img = new Image();
        img.addEventListener('error', function() {
            $s('#errfield').innerHTML = "Unable to load image " + loc;
        });
        img.addEventListener('load', function() {
            try {
                // make each canvas the image's exact size
                dims = [ img.width, img.height ];
                adjustCanvases( dims );

                // draw the image to the canvas
                ctxs[0].drawImage(img, 0, 0, img.width, img.height);

                // copy the pixels from the canvas onto rawImage
                var imageData = ctxs[0].getImageData(0, 0, dims[0], dims[1] );
                rawImage = [];
                for (var ai = 0; ai < imageData.data.length; ai+=4) {
                    // greyscale, so you only need every 4th value
                    rawImage.push(imageData.data[ai]);
                }

                // some diagnostics...
                someDiags(rawImage);
                
                $s('#errfield').innerHTML = "Loaded!";
            } catch (e) {
                $s('#errfield').innerHTML = e.message;
            }
        });
        img.crossOrigin = "anonymous";
        img.src = loc;
    }

    function generateLetter(lettervalue){
        console.log("===Step 1b: Generating an image from ", lettervalue, "===");

        // draw the letter
        dims = [256, 256] ;
        ltter = Bagel.letterImg( dims, lettervalue )
        ltter = ltter.map( x => x.magnitude() );

        // make each canvas the image's exact size
        adjustCanvases( dims );

        // draw the matrix to the canvas
        fillCanvasFrom( ctxs[0], ltter );

        // copy the pixels onto rawImage
        rawImage = ltter.slice();

        // some diagnostics...
        someDiags(rawImage);
    }

    function generatePlaid(frequencyvalues){
        console.log("===Step 1b: Generating a plaid from ", frequencyvalues, "===");
        
        // check the parameters
        try {
            var frqs = frequencyvalues.split(/[,;]/);
            if ((frqs.length<4)||(frqs.length>5)) {throw new Error("You must provide four even frequencies for 0, 45, 90 and 135 degrees, e.g., `40,0,10,0` (optional: the radius of the envelop, e.g., 40,0,10,0;70)")}
            frqs  = frqs.map( x => Math.round(parseInt(x, 10)/1)*1 );
            if (frqs.length == 5)  { var radius = frqs[4]; frqs.pop(); }
            else { var radius = (256 * 40) }; // i.e., no visible envelop
        } catch (e) {
            $s('#errfield').innerHTML = e.message;
            return;
        }

        // do the plaid
        dims = [256, 256] ;
        var angles = [Math.PI/2, 3*Math.PI/4, Math.PI, 5*Math.PI/4];
        plaid = Bagel.plaidImg( dims, angles, frqs, radius )
        plaid = plaid.map( x => x.magnitude() );

        // make each canvas the image's exact size
        adjustCanvases( dims );

        // draw the image to the canvas
        fillCanvasFrom( ctxs[0], plaid );

        // copy the pixels onto rawImage
        rawImage = plaid.slice();

        // some diagnostics...
        someDiags(rawImage);
    }


    //-----------------------------------------------------
    // Make the transformation into newImage and show in canvas1
    //-----------------------------------------------------
    function transformAction(arg, arglist, noiseprop, nselist) {
        console.log("===Step 2: Tranforming the image with args", arg, " along with ", arglist, "-- noise proportion:", noiseprop, nselist, "===");
        // convert the image to complex numbers.
        if ( rawImage[0].constructor.name != "Complex") {
            rawImage = rawImage.map( x => new Bagel.Complex(x, 0.) );
        }
        try {
            // compute newImage
            switch(arg) {
                case "0": 
                    console.log("   ...returning the raw image");
                    newImage = rawImage.slice();
                    // adding noise?
                    noiseprop = parseFloat(noiseprop, 10);
                    if (isNaN(noiseprop)||(noiseprop == 0)) {break;}
                    if ((noiseprop > 1)||(noiseprop <= 0))  {throw new Error("You must provide proportion of noise, e.g., 0.5.")}
                    var noisestd = Math.max(...newImage.map( x => x.magnitude() ) ) /4;
                    tempImage = Bagel.gaussianNoisyImg( dims, 1 ) ;
                    
                    console.log(nselist == "");
                    if ((nselist != "")) {
                        var nze = nselist.split(/[,;:]/);
                        if (nze.length<2) {throw new Error("You must provide filter type and arg(s), e.g., `f:2`")}
                        nzetype = nze.shift();
                        nzearg  = nze.map( x => parseFloat(x, 10) );
                        if (nzetype == 'f') {
                            tempFilter = Bagel.oneOverFnFilter(dims, nzearg[0] );
                            tempImage = Bagel.invFFT2D( Bagel.elementWiseTimes( Bagel.FFT2D( tempImage ), tempFilter ) );
                            tt = tempImage.map ( x=> x.magnitude() );
                            tempImage = Bagel.elementWisePlus( tempImage, -Math.min.apply(null,tt) );
                            tempImage = Bagel.elementWiseTimes( tempImage, 5/(Math.max.apply(null,tt)-Math.min.apply(null,tt))  );
                        };
                    }
                    tempImage = Bagel.elementWiseTimes( tempImage, noisestd );

                    tempImage = Bagel.elementWiseTimes( tempImage, noiseprop );
                    newImage  = Bagel.elementWiseTimes( newImage, 1-noiseprop );
                    newImage  = Bagel.elementWisePlus ( newImage, tempImage );
                    break;
                case "1": 
                    console.log("   ...returning the Fourier of image");
                    newImage = Bagel.shiftQuadrants( Bagel.FFT2D( rawImage ) );
                    break;
                
                // donut filter functions
                case "2":
                    console.log("   ...returning a band filter");
                    var bds = arglist.split(/[,;]/);
                    if (bds.length<2) {throw new Error("You must provide two radii, e.g., `3,12`, or `,12`")}
                    bds  = bds.map( x => parseInt(x, 10) );
                    newImage = Bagel.donutFilter( dims, bds[0], bds[1] );
                    break;
                case "3":
                    console.log("   ...merging a band filter (2) with the Fourier of image (1)");
                    var bds = arglist.split(/[,;]/);
                    if (bds.length<2) {throw new Error("You must provide two radii, e.g., `3,12`, or `,12`")}
                    bds = bds.map( x => parseInt(x, 10) );
                    filter   = Bagel.donutFilter( dims, bds[0], bds[1] );
                    newImage = Bagel.FFT2D( rawImage );
                    newImage = Bagel.shiftQuadrants( newImage );
                    newImage = Bagel.elementWiseTimes(filter, newImage);
                    break;
                case "4":
                    console.log("   ...untransforming the band-filtered image (3)");
                    var bds = arglist.split(/[,;]/);
                    if (bds.length<2) {throw new Error("You must provide two radii, e.g., `3,12`, or `,12`")}
                    bds = bds.map( x => parseInt(x, 10) );
                    filter   = Bagel.donutFilter( dims, bds[0], bds[1] );
                    newImage = Bagel.FFT2D( rawImage );
                    newImage = Bagel.shiftQuadrants( newImage );
                    newImage = Bagel.elementWiseTimes(filter, newImage);
                    newImage = Bagel.shiftQuadrants( newImage );
                    newImage = Bagel.invFFT2D( newImage );
                    // adding noise?
                    noiseprop = parseFloat(noiseprop, 10);
                    if (isNaN(noiseprop)||(noiseprop == 0)) {break;}
                    var noisestd = Math.max(...newImage.map( x => x.magnitude() ) ) /4;
                    if ((noiseprop > 1)||(noiseprop <= 0)) {throw new Error("You must provide proportion of noise, e.g., 0.5.")}
                    tempImage = Bagel.gaussianNoisyImg( dims, noisestd ) ;
                    tempImage = Bagel.elementWiseTimes( tempImage, noiseprop );
                    newImage  = Bagel.elementWiseTimes( newImage, 1-noiseprop );
                    newImage  = Bagel.elementWisePlus ( newImage, tempImage );
                    break;

                // gaussian filter functions
                case "5":
                    console.log("   ...returning a gauss filter");
                    var bds = arglist.split(/[,;]/);
                    if (bds == "") {throw new Error("You must provide one sigma, e.g., 12")}
                    bds  = bds.map( x => parseFloat(x, 10) );
                    newImage = Bagel.gaussFilter( dims, bds[0] );
                    break;
                case "6":
                    console.log("   ...merging a gauss filter (5) with the Fourier of image (1)");
                    var bds = arglist.split(/[,;]/);
                    if (bds == "") {throw new Error("You must provide one sigma")}
                    bds = bds.map( x => parseFloat(x, 10) );
                    filter   = Bagel.gaussFilter( dims, bds[0] );
                    newImage = Bagel.FFT2D( rawImage );
                    newImage = Bagel.shiftQuadrants( newImage );
                    newImage = Bagel.elementWiseTimes(filter, newImage);
                    break;
                case "7":
                    console.log("   ...untransforming the gauss-filtered image (6)");
                    var bds = arglist.split(/[,;]/);
                    if (bds == "") {throw new Error("You must provide one sigma")}
                    bds = bds.map( x => parseFloat(x, 10) );
                    filter   = Bagel.gaussFilter( dims, bds[0] );
                    newImage = Bagel.FFT2D( rawImage );
                    newImage = Bagel.shiftQuadrants( newImage );
                    newImage = Bagel.elementWiseTimes(filter, newImage);
                    newImage = Bagel.shiftQuadrants( newImage );
                    newImage = Bagel.invFFT2D( newImage );
                    // adding noise?
                    noiseprop = parseFloat(noiseprop, 10);
                    if (isNaN(noiseprop)||(noiseprop == 0)) {break;}
                    var noisestd = Math.max(...newImage.map( x => x.magnitude() ) ) /4;
                    if ((noiseprop > 1)||(noiseprop <= 0)) {throw new Error("You must provide proportion of noise, e.g., 0.5.")}
                    tempImage = Bagel.gaussianNoisyImg( dims, noisestd ) ;
                    tempImage = Bagel.elementWiseTimes( tempImage, noiseprop );
                    newImage  = Bagel.elementWiseTimes( newImage, 1-noiseprop );
                    newImage  = Bagel.elementWisePlus ( newImage, tempImage );
                    break;

                // bagel filter functions
                case "8":
                    console.log("   ...returning a bagel filter");
                    var dta = arglist.split(/[,;:]/);
                    if ((dta[0] != "r")&&(dta[0] != "m")) {throw new Error("You must provide a method r (random) or m (manually) as first entry")}
                    var method = dta.shift();
                    dta = dta.map( x => parseInt(x, 10) );
                    if (method == "m") {
                        var k = dta.shift();
                        var mus = dta.slice();
                        var sgs = mus.map( x => x / (2 * k) ); // heuristic
                    } else if (method = "r") {
                        var b = dta.shift();
                        var k = dta.shift();
                        var res = Bagel.bagelSampler( dims, b, k );
                        var mus = res[0];
                        var sgs = mus.map( x => x / (2 * k) ); // heuristic
                    }
                    newImage = Bagel.bagelFilter( dims, mus, sgs );
                    break;
                case "9":
                    console.log("   ...merging a bager filter (5) with the Fourier of image (1)");
                    var dta = arglist.split(/[,;:]/);
                    if ((dta[0] != "r")&&(dta[0] != "m")) {throw new Error("You must provide a method r (random) or m (manually) as first entry")}
                    var method = dta.shift();
                    dta = dta.map( x => parseInt(x, 10) );
                    if (method == "m") {
                        var k = dta.shift();
                        var mus = dta.slice();
                        var sgs = mus.map( x => x / (2 * k) ); // heuristic
                    } else if (method = "r") {
                        var b = dta.shift();
                        var k = dta.shift();
                        var res = Bagel.bagelSampler( dims, b, k );
                        var mus = res[0];
                        var sgs = mus.map( x => x / (2 * k) ); // heuristic
                    }
                    filter = Bagel.bagelFilter( dims, mus, sgs );
                    newImage = Bagel.FFT2D( rawImage );
                    newImage = Bagel.shiftQuadrants( newImage );
                    newImage = Bagel.elementWiseTimes(filter, newImage);
                    break;
                case "10":
                    console.log("   ...untransforming the bagel-filtered image (9)");
                    var dta = arglist.split(/[,;:]/);
                    if ((dta[0] != "r")&&(dta[0] != "m")) {throw new Error("You must provide a method r (random) or m (manually) as first entry")}
                    var method = dta.shift();
                    dta = dta.map( x => parseInt(x, 10) );
                    if (method == "m") {
                        var k = dta.shift();
                        var mus = dta.slice();
                        var sgs = mus.map( x => x / (2 * k) ); // heuristic
                    } else if (method = "r") {
                        var b = dta.shift();
                        var k = dta.shift();
                        var res = Bagel.bagelSampler( dims, b, k );
                        var mus = res[0];
                        var sgs = mus.map( x => x / (2 * k) ); // heuristic
                    }
                    filter = Bagel.bagelFilter( dims, mus, sgs );
                    newImage = Bagel.FFT2D( rawImage );
                    newImage = Bagel.shiftQuadrants ( newImage );
                    newImage = Bagel.elementWiseTimes(filter, newImage);
                    
                    newImage = Bagel.shiftQuadrants( newImage );
                    newImage = Bagel.invFFT2D( newImage );
                    // adding noise?
                    noiseprop = parseFloat(noiseprop, 10);
                    if (isNaN(noiseprop)||(noiseprop == 0)) {break;}
                    var noisestd = Math.max(...newImage.map( x => x.magnitude() ) ) /4;
                    if ((noiseprop > 1)||(noiseprop <= 0)) {throw new Error("You must provide proportion of noise, e.g., 0.5.")}
                    tempImage = Bagel.gaussianNoisyImg( dims, noisestd ) ;
                    tempImage = Bagel.elementWiseTimes( tempImage, noiseprop );
                    newImage  = Bagel.elementWiseTimes( newImage, 1-noiseprop );
                    newImage  = Bagel.elementWisePlus ( newImage, tempImage );
                    break;

                case "99": 
                    console.log("   ...inverting the Fourier of image (1)");
                    newImage = Bagel.FFT2D(rawImage);
                    newImage = Bagel.invFFT2D( newImage );
                    break;
                default: 
                    throw new Error("   ...You must provide argument between 0 and 10 or 99.")
                    break;
            }

            // setting if log magnitude is desired
            if ($s("#logscale").checked) { // log magnitude
                Math.squash = function(x) {return Math.log(x)};
                var shift = 1; //1 with log; 0 with identity
            } else { // magnitude as is
                Math.squash = function(x) {return x};
                var shift = 0; //1 with log; 0 with identity    
            }

            console.log("   ...normalizing the amplitudes");
            var squashedImage = newImage.map( x => Math.squash( x.magnitude() + shift) );
            // get the largest magnitude
            var maxMagnitude = Math.max(...squashedImage );
            // some diagnostics...
            someDiags(squashedImage);

            // draw the pixels
            console.log("   ...dumping the image on the screen");
            squashedImage = Bagel.elementWiseTimes( squashedImage, 255/maxMagnitude );        
            fillCanvasFrom( ctxs[1], squashedImage );

            $s('#errfield').innerHTML = "Tranformed!";
        } catch (e) {
            $s('#errfield').innerHTML = e.message;
        }
    }


    //-----------------------------------------------------
    // Make the transformation into newImage and show in canvas2
    //-----------------------------------------------------
    function contrastAction( method, range ) {
        // Rescale the range of pixel gray colors around 0.5
        // Two techniques are implemented: (a) a simple (r)escaling where
        // pixels are rescaled to cover a range from 0.5-range/2 to 0.5+range/2.
        // (b) the m(edian) rescaling whereby, as with (r), the pixels are
        // rescaled to cover a range from 0.5-range/2 to 0.5+range/2 and then
        // the pixels are shifted up or down so that the median equals 0.5.
        console.log("===Step 3: Contrasting the image with args", method);
        try {
            // compute cstImage
            tmpImage = newImage.map( x => x.magnitude() );
            minImg = Math.min.apply(null,tmpImage);
            maxImg = Math.max.apply(null,tmpImage);
            console.log( "=> Min, max  = ", [minImg, maxImg] );

            switch( method ) {
                case "r": 
                    console.log("   ...contrasting for a mean of 0.5 with a range of ", range);
                    range = parseFloat(range, 10);
                    if ((isNaN(range))||(range >1)||(range <0)) {throw new Error("You must provide gray range (between 0 and 1)")}
                    cstImage = Bagel.elementWisePlus( tmpImage, -minImg);
                    cstImage = Bagel.elementWiseTimes( cstImage, range/(maxImg-minImg) );
                    cstImage = Bagel.elementWisePlus( cstImage, 0.5-range/2 );
                    break;
                case "m": 
                    console.log("   ...contrasting for a mean of 0.5 with a range of ", range);
                    range = parseFloat(range, 10);
                    if ((isNaN(range))||(range >1)||(range <0)) {throw new Error("You must provide gray range (between 0 and 1)")}
                    cstImage = Bagel.elementWisePlus( tmpImage, -minImg);
                    cstImage = Bagel.elementWiseTimes( cstImage, range/(maxImg-minImg) );
                    cstImage = Bagel.elementWisePlus( cstImage, 0.5-range/2 );

                    mn = Math.min.apply(null,cstImage);
                    mx = Math.max.apply(null,cstImage);
                    md = Bagel.median( cstImage );
                    lg = 0.5-md;
                    console.log("Min, Max, Median, Lag:", mn, mx, md, lg);
                    cstImage = Bagel.elementWisePlus( cstImage, lg );
                    break;
                case "c": 
                    throw new Error("corner not implemented yet...");
                    break;
                default: 
                    throw new Error("   ...You must provide argument r(escale), m(edian rescale) or c(orners rescale).")
                    break;
            }

            // some diagnostics...
            someDiags( cstImage );

            // draw the pixels
            console.log("   ...dumping the image on the screen");
            cstImage = Bagel.elementWiseTimes( cstImage, 255 );        
            fillCanvasFrom( ctxs[2], cstImage );

            $s('#errfield').innerHTML = "Contrasted!";
        } catch (e) {
            $s('#errfield').innerHTML = e.message;
        }
    }



    //-----------------------------------------------------
    // managing callback functions
    //-----------------------------------------------------

    function disableButtons(callback) {
        $s('#draw-cs-btn').disabled     = true;
        $s('#draw-grace-btn').disabled  = true;
        $s('#draw-lettre-btn').disabled = true;
        $s('#draw-plaid-btn').disabled  = true;
        $s('#transform-btn').disabled   = true;
        $s('#method-btn').disabled       = true;

        setTimeout(callback, 50); // 50 ms for the UI to update
        enableButtons();
    }

    function enableButtons() {
        $s('#draw-cs-btn').disabled     = false;
        $s('#draw-grace-btn').disabled  = false;
        $s('#draw-lettre-btn').disabled = false;
        $s('#draw-plaid-btn').disabled  = false;
        $s('#transform-btn').disabled   = false;
        $s('#method-btn').disabled       = false;
    }


    //-----------------------------------------------------
    // helper functions
    //-----------------------------------------------------

    function $s(id) { // for convenience
        if (id.charAt(0) !== '#') return false;
        return document.getElementById(id.substring(1));
    }


    //-----------------------------------------------------
    // All done! The init is the only public function
    //-----------------------------------------------------
    return { 
        init: initBagelImageAnalysis 
    };
})();

window.addEventListener( 'load', BagelImageAnalysis.init );
