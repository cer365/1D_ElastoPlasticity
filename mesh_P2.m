
function [coord,elem,surf,neumann,Q]=mesh_P2(level,size_x)

% =========================================================================
%
%  This function creates triangular mesh for P2 elements
%
%  input data:
%    level     - an integer defining a density of a uniform mesh
%    size_x    - size of the body in directions x (integer)
%    body=(0,size_x)
%
%  output data:
%    coord   - coordinates of the nodes, size(coord)=(1,n_n) where n_n is a
%              number of nodes
%    elem    - array containing numbers of nodes defining each element, 
%              size(elem)=(3,n_e), n_e = number of elements
%    surf    - array containing numbers of nodes defining each surface element, 
%              size(surf)=(3,n_s), n_s = number of surface elements
%    neumann - array containing numbers of nodes defining each surface element, 
%              size(neuman)=(3,n_e_s). The surface is the following side 
%              of the body: (0,size_x), where the nonhomogeneous
%              Neumann boundary condition is considered.
%    Q       - logical array indicating the nodes where the Dirichlet
%              boundary condition is considered, size(Q)=(1,n_n)
%
% ======================================================================
%

%
% numbers of segments, nodes and elements
%

  N_x = size_x*2^level;       % number of segments in x direction
  
  % 
  n_n = 2*N_x+1;        % total number of nodes

%
% C - 1D auxilliary array that contains node numbers and that is important 
% for the mesh construction. 
%
  C=zeros(1,n_n);
  C(1,:) = 1:2*N_x+1;
  
%
% coordinates of nodes
%         
  % the required array of coordinates, size(coord)=(1,n_n)   
  coord=linspace(0,size_x,n_n);  

% 
% construction of the array elem
%
  % ordering of the nodes creating the unit cube:
  %  V1 -> [0 0], V2 -> [1 0]
  %  V1,V2 are logical 1D arrays which enable to select appropriate
  %  nodes from the array C.

  V1=false(1,n_n);
  V1(1,1:2:n_n-2)=1;
  %
  V2=false(1,n_n);
  V2(1,3:2:n_n)=1; 
  
  V12=false(1,n_n);
  V12(1,2:2:n_n-1)=1;

  % used division of a line:   
  %   V1 V2 V12        
  % the array elem, size(elem)=(3,n_e)          
  elem=[C(V1); C(V2); C(V12)];     
 
%
% Surface of the body - the array "surf"
%
  
  % the array "surf"
  surf = elem ;
  
%
% Boundary conditions
%
  
  % nonhomogeneous Neumann boundary conditions
  neumann=0;          

  % logical array indicating the nodes with the Dirichlet boundary cond.
  Q = (coord>0 & coord<size_x) ;
  
end
