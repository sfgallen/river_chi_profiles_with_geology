function [smoothEl] = smoothChannelZ(elev,win,cs)

% smooths channel data with filtfilt.
% inputs: dist = elevation in meters,  win = window size in meters,
% cs = cellsize of grid in meters

win2 = 2.*round(((win/cs)+1)/2)-1; % round to the nearest odd integer
win2 = min(win2,length(elev));

if win2*3 >= length(elev) % just do fastsmooth filter if the window length is too small.  Because the data vector needs to be >3x the length of window
    smoothEl = fastsmooth(elev,win2,3,1);
else
    smoothEl = filtfilt(ones(1,win2)/win2,1,elev);
end

smoothEl(1) = elev(1);
smoothEl(end) = elev(end);

if length(elev)>=3
    if smoothEl(end-1) < smoothEl(end)
        smoothEl(end-1) = (smoothEl(end)+smoothEl(end-2))/2;  % if point immediately upstream of confluence is below elevation of confluence, interpolate its elevation linearly between confluence and second to last point before confluence
    end
end

end