crita = 1e5;            % threshold drainage area for channel head initiation
mn = 0.45;              % m/n or concavity index
Ao = 1;               % reference drainage area for chi-analysis
smoWin = 500;           % size of window (in map units) used to smooth elevation data

addpath('C:\Users\sfgallen\Documents\topo_toolbox\topotoolbox-master');

DEM = GRIDobj('Jayuya_DEM_sm.tif');
DEM = fillsinks(DEM);
GEO = GRIDobj('Jayuya_GEO_sm.tif');

% set nan's
DEM.Z(DEM.Z <= -9999) = nan;
GEO.Z(GEO.Z <=-9999) = nan;

% make sure grids are the same size
[DEM,GEO] = largest_overlapping_extent(DEM,GEO);

cs = DEM.cellsize;

% Calculate flow object (FD)
FD  = FLOWobj(DEM,'preprocess','carve');

% Create stream network object (S)
S  = STREAMobj(FD,'minarea',crita/(cs)^2);

A = flowacc(FD).*cs^2;

[ChiGrid, chiS] = CalculateChi(S, A, Ao, mn);

geoS = GEO.Z(S.IXgrid);
geoIDs = unique(geoS(~isnan(geoS)));
cols = jet(length(geoIDs));

ordList = S.orderednanlist;
strmBreaks = find(isnan(ordList));

% plot all of the river river profile data as thin gray lines
hh = waitbar(0,'Smoothing elevation data for the stream network...');
id1 = 0;
SmoZ = double(DEM.Z(S.IXgrid));
for j = 1:length(strmBreaks)
    strmInds = ordList(id1+1:strmBreaks(j)-1);
    SmoZ(strmInds) = smoothChannelZ(SmoZ(strmInds),smoWin,cs);
    id1 = strmBreaks(j);
    ff = j/length(strmBreaks);
    waitbar(ff,hh);
end
close(hh)

figure; hold on

Sd = S.distance;
GridID = S.IXgrid;

id1 = 0;
strmNum = 0;
hh = waitbar(0,'plotting streams with geology...');
for j = 1:length(strmBreaks)
    strmInds = ordList(id1+1:strmBreaks(j)-1);
    gridInds = GridID(strmInds);
    
    sA = A.Z(gridInds);
    sD = Sd(strmInds);
    sZ = SmoZ(strmInds);
    sC = chiS(strmInds);
    sG = geoS(strmInds);
    
    inds = 1;
    
    % this is the block of code that plots the geology
    for k = 1:length(sG)-1
        if k == length(sG)-1
            colInd = find(geoIDs == sG(k));
            subplot(2,1,1)
            plot(sD(inds)./1000,sZ(inds)./1000,'Color',cols(colInd,:),'LineWidth',0.70); hold on
            subplot(2,1,2)
            plot(sC(inds),sZ(inds)./1000,'Color',cols(colInd,:),'LineWidth',0.70); hold on
            inds = 1;
        elseif sG(k) == sG(k+1)
            inds = [inds,k];
        elseif sG(k) ~= sG(k+1)
            inds = [inds,k];
            colInd = find(geoIDs == sG(k));
            subplot(2,1,1)
            plot(sD(inds)./1000,sZ(inds)./1000,'Color',cols(colInd,:),'LineWidth',0.70); hold on
            subplot(2,1,2)
            plot(sC(inds),sZ(inds)./1000,'Color',cols(colInd,:),'LineWidth',0.70); hold on
            inds = k;
        end
    end
    strmNum = strmNum + 1;
    sNumG.Z(gridInds) = strmNum;
    
    id1 = strmBreaks(j);
    ff = j/length(strmBreaks);
    waitbar(ff,hh);
end
close(hh);

subplot(2,1,1)
xlabel('distance (km)'); ylabel('elevation (km)');
subplot(2,1,2)
xlabel('\chi (m)'); ylabel('elevation (km)');