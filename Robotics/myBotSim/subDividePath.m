function subDivPath = subDividePath( path, subDiv)
%SUBDIVIDEPATH subdivides a supplied path of points into a finer
% path of points, such that the line between every two subsequent points  
% is subdivided into subDiv points.
subDivPath = path(1,:);
for i = 1:size(path,1)-1
    xvals = linspace(path(i,1),path(i+1,1),subDiv+2);
    yvals = linspace(path(i,2),path(i+1,2),subDiv+2);
    subDivPath = [subDivPath;[xvals(2:end);yvals(2:end)]'];
end
end

