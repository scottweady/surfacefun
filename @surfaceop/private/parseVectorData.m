function [data, outputType] = parseVectorData(data)
%PARSEVECTORDATA   Convert vector PDE data to scalar component columns.

outputType = 'surfacefun';

if ( isa(data, 'surfacefunv') )
    data = horzcat(data.components{:});
    outputType = 'surfacefunv';
elseif ( isnumeric(data) && isvector(data) && numel(data) >= 2 && ~isscalar(data) )
    data = reshape(data, 1, numel(data));
    outputType = 'surfacefunv';
end

end
