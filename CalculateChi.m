function [ChiGrid, ChiStreams] = CalculateChi(S, A, Ao, mn)
%
% CalculateChi.m calculates chi using a TopoToolBox STREAMobj (S), Flow
% accumulation GRIDobj (A), reference drainage area (Ao) and the m/n ratio
% (mn).
%
% The outputs are a GRIDobj of chi (ChiGrid) and chi for the topologically
% ordered stream network (ChiStreams)
%
% Author: Sean F. Gallen
% Date Modified: 12/21/2015

% get variables ready for chi integration
ChiStreams = zeros(size(S.distance));
Six = S.ix;
Sixc = S.ixc;
Sx = S.distance;
Sa = (Ao./(A.Z(S.IXgrid))).^mn;

h = waitbar(0,'calculating \chi for the stream network...');
% calculating chi for the entire river network base on netcumtrapz from
% topotoolbox
for lp = numel(Six):-1:1;
    ChiStreams(Six(lp)) = ChiStreams(Sixc(lp)) + (Sa(Sixc(lp))+(Sa(Six(lp))-Sa(Sixc(lp)))/2) *(abs(Sx(Sixc(lp))-Sx(Six(lp))));
    f = (numel(Six)+1 - lp)/numel(Six);
    waitbar(f,h);
end
close(h);

ChiGrid = A;
ChiGrid.Z = nan(size(A.Z));
ChiGrid.Z(S.IXgrid) = ChiStreams;
end