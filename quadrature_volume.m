function [Xi, WF] = quadrature_volume(elem_type)

%--------------------------------------------------------------------------
% This function specifies a numerical quadrature for wire integration,
% depending on a chosen finite element. The quadratures suggested
% below can be simply replaced by another ones.
%
%  input data: 
%       elem_type - the type of Lagrange finite elements
%
%  output data:
%       Xi - local coordinates of quadrature points, size(Xi)=(1,n_q)
%       WF - weight factors, size(WF)=(1,n_q_s)
%--------------------------------------------------------------------------

switch(elem_type)
    
  case {'P1'}      
    % - surface is created by the reference line with coordinates:
    %               -1, 1
    % - 1-point quadrature rule, i.e., n_q=1
    Xi=0; WF=2; 
   
  case {'P2'}
    % - surface is created by the reference line with coordinates:
    %               -1, 1
    % - 2-point quadrature rule, i.e., n_q=2
    pt = 1/sqrt(3);
    Xi=[-pt,pt];
    WF=[1,1]; 
      
  otherwise; disp('Bad choise of element type');
end
