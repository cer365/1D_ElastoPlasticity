function draw_quantity(coord,U,Q_node)

% =========================================================================
%
%  This function depicts prescribed nodal quantity
%
%  input data:
%    coord     - coordinates of the nodes, size(coord)=(1,n_n) where n_n is a
%                number of nodes
%    U         - nodal displacements, size(U)=(1,n_n) to catch deformed shape
%                if the deformed shape is not required then set 0*U
%    Q_node    - prescribed nodal quantity, size(Q_node)=(1,n_n)
%    elem_type - the type of finite elements; available choices:
%                'P1', 'P2'
%    size_x   - size of the body in x direction (integer)
%                size_hole < size_xy
%    body=(0,size_x)
%
% ======================================================================
%

  % preprocesing for drawing
  width_line = 0.02;
  n_n = size(coord,2);
  coord_draw = [coord, coord; -(width_line/2)*ones(1,n_n), (width_line/2)*ones(1,n_n)];
  aux = [1:n_n; n_n+1:2*n_n];
  elem_draw = [aux(1,1:n_n-1); aux(1,2:n_n); aux(2,2:n_n); aux(2,1:n_n-1)];
  U_draw = [zeros(1,2*n_n); U, U];
  Q_node_draw = [Q_node, Q_node];

  figure;
  hold on;
  
  % visualization of the quantity

  s = patch('Faces',elem_draw(1:4,:)','Vertices',coord_draw'+U_draw',...
    'FaceVertexCData',Q_node_draw','FaceColor','interp','EdgeColor','none'); 

  colorbar;

  %
  box on
  view(2);      % standard view ve 2D
  axis equal;   % real ratios
  hold off;
  %axis off;
  axis([0 1 -0.2 0.2]);
end