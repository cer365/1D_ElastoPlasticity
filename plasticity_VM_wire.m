% =========================================================================
%
%  This program triggers an assembly test for a 1D elastoplastic body. It
%  is considered the von Mises yield criterion and a linear isotropic
%  hardening. The tangent stiffness matrices
%  are computed in each time step. One can set optionally 2 types of finite
%  elements, levels of mesh density and many other parameters.
%
% ======================================================================
%

%
% The main input data 
%

  elem_type='P1'; % the type of finite elements; available choices:
                  % 'P1', 'P2'
  level=10;        % a nonnegative integer defining mesh density
  
  % values of elastic material parameters
  Ec= 206900;                        % Young's modulus 
  
  % values of plastic material parematers
  H_m=10000;
  sigma_Y=450*sqrt(2/3);
  
  % constant tranction on the top of the body in each direction
  traction_force = 1000;
           
%
% Mesh generation
%

  % geometrical parameters (choose only integers)
  size_x = 1;      % size of the wire in direction x 
 
                     
  % the mesh generation depending prescribed finite elements        
  switch(elem_type)
    case 'P1'
        [COORD,ELEM,SURF,NEUMANN,Q]=mesh_P1(level,size_x);
        fprintf('P1 elements: \n')
    case 'P2'
        [COORD,ELEM,SURF,NEUMANN,Q]=mesh_P2(level,size_x);
        fprintf('P2 elements: \n')
    otherwise
          disp('bad choice of element type');
  end  

%
% Data from the reference element
%
  
  % quadrature points and weights for volume integration
  [Xi, WF]     = quadrature_volume(elem_type);
  
  % local basis functions and their derivatives for volume 
  [HatP,DHatP1] = local_basis_volume(elem_type, Xi);

%
% Number of nodes, elements and integration points + print
%
  n_n=size(COORD,2);          % number of nodes
  n_unknown=length(COORD(Q)); % number of unknowns
  n_e=size(ELEM,2);           % number of elements
  n_q=length(WF);             % number of quadratic points
  n_int = n_e*n_q ;           % total number of integrations points 
  % 
  fprintf('number of nodes =%d ',n_n);  
  fprintf('\n');   
  fprintf('number of unknowns =%d ',n_unknown);
  fprintf('\n');   
  fprintf('number of elements =%d ',n_e);
  fprintf('\n');   
  fprintf('number of integration points =%d ',n_int);
  fprintf('\n');   

%
% Assembling of the elastic stiffness matrix
%
 
  % values of elastic material parameters at integration points
  Ec =Ec*ones(1,n_int);
  
  % stiffness matrix assembly and the assembly time
  tic;     
  [K_elast,B,WEIGHT,iD,jD,D_elast]=elastic_stiffness_matrix(ELEM,COORD,...
                            Ec,DHatP1,WF);  
  assembly_elast_time=toc; 
  fprintf('step =%d ',1);
  fprintf('\n');   
  fprintf('  time spent on K_elast:  %6.1e seconds, ',assembly_elast_time);
  fprintf('\n');   

%
% Assembling of the vector of traction forces
%  
  
  % values of the density function f_t at surface integration points
  n_e_s=size(ELEM,2);    % number of surface elements
  n_q_s=length(WF);      % number of quadrature points on a surface element
  n_int_s=n_e_s*n_q_s ;  % number of integration points on the surface
                         % (on the upper side of the body)
  f_t_int=traction_force'*ones(1,n_int_s); % size(f_V_int)=(1,n_int_s)
  
  % assembling of the vector of traction (surface) forces
  f_t=vector_traction(ELEM,COORD,f_t_int,HatP,DHatP1,WF);
  
%
% Loading process and Newton's solver
%

  % number of load steps and values of load factors
  d_zeta=1/10;              % constant load increment
  zeta=[0:d_zeta:1, (1-d_zeta):-d_zeta:(-1.2), (-1.2+d_zeta):d_zeta:0];
                            % sequence of load factors
  n_step = length(zeta);    % number of load steps
  alpha=zeros(1,n_step);    % work of external forces  
  
  % values of plastic material parematers at integration points
  H_m=H_m*ones(1,n_int);
  sigma_Y=sigma_Y*ones(1,n_int);
  
  % initialization
  U = zeros(1,n_n);
  U_zeros = zeros(1,n_n);
  dU = zeros(1,n_n) ;        % Newton's increment (in displacement)
  U_old = zeros(1,n_n) ; 
  F = zeros(1,n_n) ;         % vector of internal forces
  E = zeros(1,n_int);        % strain tensors at integration points  
  S_old = zeros(1,n_int);    % previous stress tensors at integration points  
  S_new = zeros(1,n_int);    % previous stress tensors at integration points
  Hard_old = zeros(1,n_int); % hardening tensors at integration points
  Hard_new = zeros(1,n_int); % hardening tensors at integration points
  
  % storage of assembly time in dependence on plastic integration points
  assembly=zeros(20*n_step,2);
  assembly_step=0;   
  
  % for loop through load steps 
  for i=2:n_step            
      
    fprintf('step =%d ',i);
    fprintf('\n');     
     
    f=(zeta(i)-zeta(i-1))*f_t;     % the load vector at step i
    
    % values from previous time step
    S_old = S_new;
    Hard_old = Hard_new;
    
    % initial displacements
    %U_it=U;
    U_it = U_zeros;
      
    % Newton's solver (the semismooth Newton method)
    it=0;              % iteration number
    while true       
        
      % consistent tangent stiffness matrix 
      tic
      % strain at integration points
      E(:) = B*U_it(:) ;
      % solution of the constitutive problem
      [S,DS,IND_p]=constitutive_problem(E,S_old,Hard_old,Ec,H_m,sigma_Y);
      vD = repmat(WEIGHT,1,1).*DS ;
      D_p = sparse( iD(:),jD(:),vD(:), 1*n_int,1*n_int ) ;   
      K_tangent = K_elast+B'*(D_p-D_elast)*B;   
      assembly_time=toc;
      
      % measuring assembly dependance on plastic integration points
      n_plast=length(WEIGHT(IND_p));
      assembly_step=assembly_step+1;
      assembly(assembly_step,:)=[n_plast assembly_time];      
      
      fprintf('  time spent on K_tangent: %6.1e seconds, ',assembly_time);
      fprintf('  plastic integration points: %d (of %d), ',n_plast,n_int); 
 
      % vector of internal forces
      F(:) = B'*reshape(repmat(WEIGHT,1,1).*S(1:1,:), 1*n_int,1) ;
      
      % Newton's increment
      dU(Q) = (K_tangent(Q,Q)\(f(Q)-F(Q))')'; 
             
      % next iteration
      U_new= U_it + dU ;

      % stopping criterion 
      q1 = sqrt( dU(:)'*K_elast*dU(:) ) ;
      q2 = sqrt(  U_it(:)'*K_elast*U_it(:)  ) ;
      q3 = sqrt( U_new(:)'*K_elast*U_new(:) ) ;
      criterion = q1/(q2+q3);
      
      fprintf('  stopping criterion=%6.1e  ',criterion); 
      fprintf('\n');  
      
      % update of unknown arrays
      U_it=U_new;                                                     
            
      % test on the stopping criterion
      if  criterion < 1e-12
          break
      end
      
      % test on number of iteration
      it=it+1; 
      if  it > 50
          error('The Newton solver does not converge.')
      end         
    end%  true
     
    U = U_old + U_it;
    E(:) = B*U_it(:) ;
    [S,DS,IND_p,Hard]=constitutive_problem(E,S_old,Hard_old,Ec,H_m,sigma_Y); 
    S_new = S_old + S;
    U_old=U;
    Hard_new = Hard_old + Hard;   
    alpha(i)=f_t(Q)*U(Q)';
    
    if (i==11)|(i==21)|(i==33)|(i==45)
               
        % hardening
        Hard_node = transformation(Hard_new,ELEM,WEIGHT); 
        draw_quantity(COORD,10*U,zeros(size(Hard_node))+Hard_node); 
        colorbar off; colorbar('location','south')
        
        % stress
        Stress_node = transformation(S_new,ELEM,WEIGHT);
        draw_quantity(COORD,10*U,zeros(size(Stress_node))+Stress_node); 
        colorbar off; colorbar('location','south')
        
    end
         
  end %for

%  
% Postprocessing - visualization of selected results
%
      
  % mesh visualization
  draw_mesh(COORD)     
  
  % load path
  figure; hold on; 
  plot(alpha,zeta,'x-');
  plot(alpha(11),zeta(11),'ro'); 
  axis([-1 1.5 -1.5 1.5]);
  legend('all time steps', 'visualized time step');
  hold off
     
  % total displacements + deformed shape
  draw_quantity(COORD,10*U,U)
