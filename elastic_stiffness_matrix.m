
function  [K,B,WEIGHT,iD,jD,D]=elastic_stiffness_matrix(ELEM,COORD,...
                                            E,DHatP1,WF)
 
% =========================================================================
%
%  This function assembles the elastic stiffness matrix and some other
%  auxilliary arrays which are important for assembling of the tangent
%  stiffness matrix related to the elastoplacitic body
%
%  input data:
%    ELEM   - array containing numbers of nodes defining each element,
%             size(ELEM)=(n_p,n_e), n_e = number of elements
%    COORD  - coordinates of the nodes, size(coord)=(1,n_n) where n_n is a
%             number of nodes
%    E      - E moduli at integration points, size(shear)=(1,n_int)
%    DHatP1 - derivatives of basis functions at the quadrature points 
%             in the direction xi_1, size(HatP)=(n_p,n_q)
%    WF     - weight factors on the reference element, size(WF)=(1,n_q)
%
%  output data:
%    K      - the elastic stiffness matrix, size(K)=(2*n_n, 2*n_n)
%    B      - the strain-displacement matrix, size(B)=(3*n_int,2*n_n)
%    WEIGHT - the weight coefficients for each quadrature point, 
%             size(WEIGHT)=(1,n_int)
%
% ======================================================================
%

%
% auxilliary notation
%

  n_n=size(COORD,2);    % number of nodes including midpoints
  n_e=size(ELEM,2);     % number of elements
  n_p=size(ELEM,1);     % number of vertices per element
  n_q=length(WF);       % number of quadrature points
  n_int = n_e*n_q ;     % total number of integrations points
     
%
% Jacobian, its determinant and inverse, 
% derivative of local basis functions
% 

  % extension of the input arrays DHatP1 by replication
  % size(DHatPhi1)=(n_p,n_int)
  DHatPhi1=repmat(DHatP1,1,n_e);
  
  % coordinates of nodes defining each element
  % size(COORDe1)=(n_p,n_e)
  COORDe1=reshape(COORD(1,ELEM(:)),n_p,n_e);
  
  % coordinates of nodes around each integration point
  % size(COORDint1)=(n_p,n_int)
  COORDint1=kron(COORDe1,ones(1,n_q)); 
  
  % components of the Jacobians: size=(1,n_int)
  J11=sum(COORDint1.*DHatPhi1);
  
  % determinant of the Jacobian: size=(1,n_int)
  DET=J11;
  
  % components of the inverse to the Jacobian: size=(1,n_int)
  Jinv11 = 1./DET;
  
  % derivatives of local basis functions w.r.t the coordinates x_1:
  % size(DPhi1)=(n_p,n_int)
  DPhi1 = repmat(Jinv11,n_p,1).*DHatPhi1;
  
%
% assembling of the strain-displacement matrix
% size(B)=(1*n_int,1*n_n)
%   
  
  % values of the strain-displacement matrix B
  n_b = 1*n_p ;
  vB=zeros(n_b,n_int);        % size(vB)=(6*n_p,n_int)
  vB(1:n_b,:)=DPhi1; 

  % i-th and j-th indices of B: size(iB)=size(jB)=(1*n_p,n_int)
  AUX=reshape(1:1*n_int,1,n_int);
  iB=repmat(AUX,1*n_p,1);  

  AUX1 = 1*(1:n_p);     % size(AUX1)=(1,n_p)
  AUX2 = 0*ones(1,n_p); % size(AUX2)=(1,n_p)
  AUX3 = ELEM((AUX1(:))',:)-kron(ones(1,n_e),AUX2(:));
                              % size(AUX3)=(2*n_p,n_p)
  jB=kron(AUX3,ones(1,n_q));
  
  % the sparse strain-displacement matrix B
  B=sparse(iB(:),jB(:),vB(:), 1*n_int,1*n_n);  

%  
% assembling of the elastic stress-strain matrix 
% size(D)=(1*n_int,1*n_int)
%
  
  % elastic tensor at each integration point: 
  Elast=E; % size(Elast)=(1,n_int)
 
  % weight coefficients: size(WEIGHT)=(1,n_int)
  WEIGHT = abs(DET).*repmat(WF,1,n_e);
  
  % assemblinng of the sparse matrix D
  iD=repmat(AUX,1,1); 
  jD=kron(AUX,ones(1,1));
  vD=Elast.*(ones(1,1)*WEIGHT);
  D=sparse(iD,jD,vD);

%
% elastic stiffness matrix: size(K)=(1*n_n,1*n_n)
%
  K = B'*D*B ;    
 
end  % end of function