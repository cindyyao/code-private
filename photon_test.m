mglVisualAngleCoordinates(57,[16 12]);
mglClearScreen(0.5);
grating = mglMakeGrating(16,12,1.5,45,0);
% gaussian = mglMakeGaussian(16,12,1,1); 
% gabor = 255*(grating.*gaussian+1)/2;
% tex = mglCreateTexture(gabor);
gabor = 255*(grating+1)/2;
tex = mglCreateTexture(gabor);
mglBltTexture(tex,[0 0]);
% mglScreenCoordinates
% mglQuad([0; 0; 100; 100], [0; 600; 600; 0], [.5; .5; .5], 0);
% mglQuad([700; 700; 800; 800], [0; 600; 600; 0], [.5; .5; .5], 0);
% mglQuad([0; 0; 800; 800], [0; 100; 100; 0], [.5; .5; .5], 0);
% mglQuad([0; 0; 800; 800], [500; 600; 600; 500], [.5; .5; .5], 0);

mglFlush;

mglScreenCoordinates
mglClearScreen(0.5);
grating = mglMakeGrating(600,400,75,45,0);
% gaussian = mglMakeGaussian(16,12,1,1); 
% gabor = 255*(grating.*gaussian+1)/2;
% tex = mglCreateTexture(gabor);
gabor = 255*(grating+1)/2;
tex = mglCreateTexture(gabor);
mglBltTexture(tex, [400 300]);
mglFlush;


numCycles = 10; texWidthPixels = 400; obj.color = [1 1 1]; obj.backgrndcolor = [0 0 0];


grating = 255*(sign(sin(0:numCycles*2*pi/(texWidthPixels-1):numCycles*2*pi))+1)/2;
colored_grating = cat(3, ( (grating .* obj.color(1)) + round(255 .* obj.backgrndcolor(1)) ), ( (grating .* obj.color(2)) + round(255 .* obj.backgrndcolor(2)) ), ( (grating .* obj.color(3)) + round(255 .* obj.backgrndcolor(3)) ));
tex1dsquare = mglCreateTexture(colored_grating);
mglBltTexture( tex1dsquare, [400 300 nan texWidthPixels], 0, 0, 45 ); 

mglFlush();
