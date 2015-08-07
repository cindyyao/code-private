function mask = make_circular_mask(center, radius, width, height, scale)

if nargin < 5
    scale = 1;
end

cx = center(1)*scale;
cy = center(2)*scale;
radius = radius*scale;
width = width*scale;
height = height*scale;


[x,y] = meshgrid(-(cx-1):(width-cx),-(cy-1):(height-cy));
mask = ((x.^2+y.^2)<=radius^2);
end