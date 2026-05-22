function out = isempty(f)
%ISEMPTY   Test for empty SURFACEFUNV.
%   ISEMPTY(F) returns 1 if F is an empty SURFACEFUNV and 0 otherwise.

out = isempty(f.components) || all(cellfun(@isempty, f.components));

end
