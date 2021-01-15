function [GRID1, GRID2] = largest_overlapping_extent(GRID1,GRID2)
%
% largest_overlapping_extent.m take two TopoToolbox GRIDobjects and resizes
% them to their largests overlapping extent.

%% set nan values for common nan numbers
GRID1.Z(GRID1.Z<=-9999) = nan;
GRID1.Z(GRID1.Z>=9999) = nan;

% set nan values for common nan numbers
GRID2.Z(GRID1.Z<=-9999) = nan;
GRID2.Z(GRID1.Z>=9999) = nan;

%% resize grids
% make sure grids are the same size
[xMin(1,1),xMax(1,1),yMin(1,1),yMax(1,1)] = findCorners(GRID1);
[xMin(2,1),xMax(2,1),yMin(2,1),yMax(2,1)] = findCorners(GRID2);

% resizing grids to largest overlapping area
xMinP = max(xMin);
xMaxP = min(xMax);
yMinP = max(yMin);
yMaxP = min(yMax);

GRID1 = gridReSize(GRID1,xMinP,xMaxP,yMinP,yMaxP);
GRID2 = gridReSize(GRID2,xMinP,xMaxP,yMinP,yMaxP);

%% revise GRIDobj and the refmat as needed
[Ny1,Nx1] = size(GRID1.Z);
[Ny2,Nx2] = size(GRID2.Z);

% fix rows if they aren't the same length
if length(Ny1) < length(Ny2)
    l_dif = length(Ny2) - length(Ny1);
    GRID2.Z = GRID2.Z(1:end-l_dif,:);
elseif length(Ny1) > length(Ny2)
    l_dif = length(Ny1) - length(Ny2);
    GRID1.Z = GRID1.Z(1:end-l_dif,:);
end

% fix columns if they aren't the same length
if length(Nx1) < length(Nx2)
    l_dif = length(Nx2) - length(Nx1);
    GRID2.Z = GRID2.Z(:,1:end-l_dif);
elseif length(Nx1) > length(Nx2)
    l_dif = length(Nx1) - length(Nx2);
    GRID1.Z = GRID1.Z(:,1:end-l_dif);
end

GRID1.size = size(GRID1.Z);
GRID2.size = size(GRID2.Z);

% % fix refmat
% [Ny1,Nx1] = size(GRID1.Z);
% 
% GRID1.georef.SpatialRef.XWorldLimits = [0, Nx1*GRID1.cellsize];
% GRID1.georef.SpatialRef.YWorldLimits = [0, Ny1*GRID1.cellsize];

GRID2.refmat = GRID1.refmat;

end

function [xMinP,xMaxP,yMinP,yMaxP] = findCorners(inGrid)
% the findCorners function takes an acsii grid input made with the function
% 'makeGrid' and finds the UTM coorinates of the corners of the grid.
% Created: 11/06/2013
% Author: Sean F. Gallen

ncols = inGrid.size(2);
nrows = inGrid.size(1);

[yMinP,xMinP] = pix2latlon(inGrid.refmat,nrows,1);
[yMaxP,xMaxP] = pix2latlon(inGrid.refmat,1,ncols);

end

function outGrid = gridReSize(inGrid,xMinP,xMaxP,yMinP,yMaxP)
% gridReSize takes five inputs, an grid and maximum and minimum
% corrdinates you want to strink the grid to. The output is the resized
% grid.
% Created: 11/06/2013
% Author: Sean F. Gallen

ncols = inGrid.size(2);
nrows = inGrid.size(1);
R = inGrid.refmat;
[RSmaxR, RSminC] = UTMlatlon2pix(R,yMinP, xMinP);
[RSminR, RSmaxC] = UTMlatlon2pix(R,yMaxP, xMaxP);
GminR = 1;
GminC = 1;
GmaxC = ncols;
GmaxR = nrows;

%% This block of code resizes the grids according to smallest dimentions 
% and creates headers for the ascii file associated witht the resized grids

if RSminR >= GminR;
    minR = RSminR;
else
    minR = GminR;
end

if RSmaxR <= GmaxR;
    maxR = RSmaxR;
else
    maxR = GmaxR;
end

if RSminC >= GminC;
    minC = RSminC;
else
    minC = GminC;
end

if RSmaxC <= GmaxC;
    maxC = RSmaxC;
else
    maxC = GmaxC;
end
 outGrid = inGrid;
 
 minR = round(minR);
 maxR = round(maxR);
 minC = round(minC);
 maxC = round(maxC);
 
 rows = minR:maxR;
 cols = minC:maxC;
 outGrid.Z = inGrid.Z(rows,cols);
 outGrid.size = size(outGrid.Z);
 [yll,xll] = pix2UTMlatlon(inGrid.refmat,minR,minC);
 outGrid.refmat = makerefmat(xll, yll, inGrid.cellsize, -inGrid.cellsize);
 
end

function [lat, lon] = pix2UTMlatlon(R, row, col)

dx = R(2,1);
dy = R(1,2);

xll = R(3,1);
yll = R(3,2);

lat = row*dy+yll;
lon = col*dx+xll;
end

function [row, col] = UTMlatlon2pix(R, lat, lon)

dx = R(2,1);
dy = R(1,2);

xll = R(3,1);
yll = R(3,2);

row = (lat-yll)/dy;
col = (lon-xll)/dx;
end