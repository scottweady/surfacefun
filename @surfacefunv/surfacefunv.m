classdef surfacefunv
%SURFACEFUNV   A class for representing functions on surfaces.
%   SURFACEFUNV(DOM) constructs an SURFACEFUNV over the SURFACEMESH DOM.
%
%   SURFACEFUNV(F, G, H) constructs a SURFACEFUNV representing the vector
%   field (F, G, H), where F, G, and H are SURFACEFUNs with the same
%   domain.
%
%   SURFACEFUNV(F, G, H, DOM) constructs a SURFACEFUNV representing the
%   vector field (F, G, H) over the SURFACEMESH DOM, where F, G, and H are
%   function handles of the form @(x,y,z) ... or cell arrays of function
%   values at tensor-product Chebyshev nodes.

    properties

        components
        isTransposed

    end

    methods

        function f = surfacefunv(varargin)

            if ( nargin == 0 )
                % Empty SURFACEFUNV:
                f.components = {surfacefun, surfacefun, surfacefun};
            elseif ( nargin == 1 && isa(varargin{1}, 'surfacemesh') )
                % Call is: SURFACEFUNV(DOM)
                dom = varargin{1};
                f.components = {surfacefun(dom), ...
                                surfacefun(dom), ...
                                surfacefun(dom)};
            elseif ( nargin >= 2 && all(cellfun(@(v) isa(v, 'surfacefun'), varargin)) )
                % Call is: SURFACEFUNV(F1, F2, ..., FN) with N >= 2 surfacefuns
                f.components = varargin;
            elseif ( nargin >= 3 && isa(varargin{end}, 'surfacemesh') && ...
                     all(cellfun(@(v) ~isa(v, 'surfacemesh'), varargin(1:end-1))) )
                % Call is: SURFACEFUNV(F1, ..., FN, DOM)
                dom = varargin{end};
                funs = varargin(1:end-1);
                f.components = cellfun(@(g) surfacefun(g, dom), funs, ...
                                       'UniformOutput', false);
            else
                error('SURFACEFUNV:surfacefunv:invalid', ...
                        'Invalid call to surfacefunv constructor.');
            end

            f.isTransposed = false;

        end

    end

end
