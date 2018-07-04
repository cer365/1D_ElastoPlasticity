
function [S,DS,IND_p,Hard]=constitutive_problem(E,S_prev,Hard_prev,...
                                                      Ec,H_m,Y)

% =========================================================================
%
% The aim of this function is to construct constitutive and consistent 
% tangent operators at integration points 1,2,...,n_int. These operators
% are related to elastoplastic model containing the von Mises yield
% criterion and a linear isotropic hardening.
%
% Input data:
%  E         - current strain tensor, size(E)=(1,n_int)
%  S_prev    - stress tensor from the previous time step,
%              size(S_prev)=(1,n_int)
%  Hard_prev - isotropic hardening (tensor) from the previous time step
%              size(Hard_prev)=(1,n_int)
%  Ec        - material constant at integration points, size(Ec)=(1,n_int)
%  H_m       - hardening parameters at integration points, size(a)=(1,n_int)
%  Y         - initial yield stress at integration points, size(Y)=(1,n_int)
%
% Output data:
%  S     - stress tensors at integration points, size(S)=(1, n_int)
%  DS    - consistent tangent operators at integr. points,
%          size(DS)=(1,n_plast)
%  IND_p - logical array indicating integration points with plastic response, 
%          size(IND_p)=(1,n_int), 
%          n_plast=number of the points with plastic response
%  Hard  - tensors of kinematic hardening, size(Hard)=(4,n_int)
%
% =========================================================================
 
  hard_part = Y + H_m.*Hard_prev;
  S_el = Ec.*E;
  S_tr = S_prev + S_el;
  DS_tr= (Ec.*H_m)./(Ec + H_m);
  tau_eps = 1e-6;
  
%
% Evaluation of the yield criterion and specification of integration points
% with plastic response
%

  CRIT1=S_tr+hard_part;        % size(CRIT1)=(1,n_int)
  IND_p1=CRIT1<-tau_eps;       % logical array, size(IND_p1)=(1,n_int)   
  CRIT2=S_tr-hard_part;        % size(CRIT2)=(1,n_int)
  IND_p2=CRIT2> tau_eps;       % logical array, size(IND_p2)=(1,n_int)

%
% The elastic prediction of unknowns
%

  S=S_el; DS=Ec;

%
% The plastic correction at the selected integration points
%
  denom = Ec + H_m;
  Phi1 = -S_tr - hard_part;
  Phi2 =  S_tr - hard_part;
  
  S(IND_p1) = S_el(IND_p1) + Ec(IND_p1)./denom(IND_p1).*Phi1(IND_p1);
  S(IND_p2) = S_el(IND_p2) - Ec(IND_p2)./denom(IND_p2).*Phi2(IND_p2);
  
  DS(IND_p1) = DS_tr(IND_p1);
  DS(IND_p2) = DS_tr(IND_p2);
  
  IND_p = IND_p1 | IND_p2;

%     
% Update of the plastic strain and the hardening variable (optional)
%

  if nargout>3    
     Hard = zeros(size(Hard_prev));
     Hard(IND_p1) =  (ones(nnz(IND_p1),1)')./denom(IND_p1).*Phi1(IND_p1);
     Hard(IND_p2) =  (ones(nnz(IND_p2),1)')./denom(IND_p2).*Phi2(IND_p2);
  end

 end
