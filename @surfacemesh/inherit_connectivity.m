
function [x, y] = inherit_connectivity(x, y, dom)
%INHERIT_CONNECTIVITY   Inherit connectivity from a reference domain.

  edges = dom.connectivity.edges; %which nodes are connected
  edge2elem = dom.connectivity.edge2elem; %which elements an edge touches

  for ne = 1 : length(edges)

    % Find which elements 
    k1k2 = edge2elem(ne, :);
    
    % First element
    k1 = k1k2{1}(1);

    % Second element
    k2 = k1k2{1}(2);

    % Check if are connected
    if k2 ~= -1

      % Get coordinates
      x1 = x{k1}; y1 = y{k1};
      x2 = x{k2}; y2 = y{k2};

      % Edge constructor
      edges = @(x, y) {[x(1,:); y(1, :)], [x(end,:); y(end, :)], [x(:,1) y(:,1)]', [x(:,end),y(:,end)]'};
      e1 = edges(x1, y1);
      e2 = edges(x2, y2);

      % Find matching edge
      err_min = inf;
      err_flip_min = inf;

      for i = 1 : 4
        for j = 1 : 4
          
            ei = e1{i};
            ej = e2{j};

            err = max(abs(ei - ej), [], 'all');
            err_flip = max(abs(ei - fliplr(ej)), [], 'all');

            if err < err_min
              err_min = err;
              i_min = [i j];
            end

            if err_flip < err_flip_min
              err_flip_min = err_flip;
              i_flip_min = [i j];
            end

        end
      end

      if err_flip_min < err_min
          is_reversed = true;
          i = i_flip_min(1);
          j = i_flip_min(2);
      else
        is_reversed = false;
        i = i_min(1);
        j = i_min(2);
      end
         
      ei = e1{i};

      if is_reversed
        ei = fliplr(ei);
      end

      switch j
        case 1, x2(1,:) = ei(1, :);
                y2(1,:) = ei(2, :);
        case 2, x2(end,:) = ei(1, :);
                y2(end,:) = ei(2, :);
        case 3, x2(:,1) = ei(1, :)';
                y2(:,1) = ei(2, :)';
        case 4, x2(:,end) = ei(1, :)';
                y2(:,end) = ei(2, :)';
      end

      % Update
      x{k1} = x1; x{k2} = x2;
      y{k1} = y1; y{k2} = y2;

    end

  end
