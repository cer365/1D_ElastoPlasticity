function [HatP,DHatP1] = local_basis_volume(elem_type, Xi)

%--------------------------------------------------------------------------
% This function evaluates local basis functions on a surface and their 
% derivatives at prescribed quadrature points depending on a chosen
% finite elements.
%
%  input data: 
%       elem_type - the type of Lagrange finite elements
%       Xi        - coordinates of the quadrature points, size(Xi)=(1,n_q)
%
%  output data:
%       HatP   - values of basis functions at the quadrature points,
%                  size(HatP_s)=(n_p,n_q)
%       DHatP1 - derivatives of basis functions at the quadrature points 
%                  in the direction xi_1, size(DHatP1_s)=(n_p,n_q)
%       n_p    - number of basis functions on the element
%       n_q    - number of integration points within the elements
%--------------------------------------------------------------------------

xi = Xi(1,:); 

switch(elem_type)
  case 'P1'
    % - the reference line with coordinates:
    %               -1, 1
    % - n_p=2, n_q=length(xi_1)
    HatP = (1/2)*[1-xi; 1+xi] ;
    DHatP1 = [-1/2; 1/2];
    
  case 'P2'
    % - the reference line with coordinates:
    %               -1, 1
    % - coordinates of midpoint:
    %             0
    % - n_p=3, n_q=length(xi_1)
    HatP = [xi.*(xi-1)/2; xi.*(xi+1)/2; (xi+1).*(1-xi)] ;
    DHatP1 = [ xi-1/2; xi+1/2; -2*xi];
       
  otherwise; disp('Bad choise of element type');
end
