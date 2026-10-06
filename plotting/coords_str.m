function s = coords_str(xvals, yvals)
% COORDS_STR
%   Formats (x,y) pairs as pgfplots "coordinates{ ... }" content, using
%   the same "%e, %e" style as the reference document.
    parts = arrayfun(@(x,y) sprintf('(%e, %e)', x, y), ...
        xvals(:), yvals(:), 'UniformOutput', false);
    s = strjoin(parts, ' ');
end