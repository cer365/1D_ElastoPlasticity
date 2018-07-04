function draw_mesh(coord)

% =========================================================================
%
%  This function draws mesh and nodal point on the body
%
%  input data:
%    coord - coordinates of the nodes, size(coord)=(1,n_n) where n_n is a
%            number of nodes
% ======================================================================
%

  figure
  hold on
  plot( coord(1,:),zeros(size(coord)), '-b.', 'MarkerSize',10);
  axis equal;  % real ratios
  view(2);     % standard view in 2D
  hold off;
  axis off;
end