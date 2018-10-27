% CEE 530 Final Project
% Author: Haitong Wang
% 2D time-dependent FEM Diffusion problem

clear all;

% Inputs
height = 0.34;% cm; depth of the tissue layer
length = 5; % cm; length of the tissue layer
D_T = 5E-08; % cm2/sec, diffusion coefficient for the tissue
k_B = 0.122 / 3600; % sec-1, damping ratio for blood
c0 = 0.1; % g/ml, initial concentration
ct0 = 0; % initial condition for t = 0;
dt = 1; % time step size

% generate the 2d mesh tiangular mesh
% x: the x coordinate, y: y coordinate, node: element node number,
[x,y,node,numele, numnod,ndivl,ndivw] = mesh2d_diff(length, height);

% a =1;
% for i = 1: numele 
%     for j = 1:3
%     index = node(j,i);
%     xplot(a) = x(index);
%     yplot(a) = y(index);
%     a = a+1;
%     end 
% end 

% plot the mesh node 
% figure 
% plot(xplot,yplot,'-o')
% title('3 node triangular mesh')
% xlabel('length(cm)')
% ylabel('height(cm)')

% set the two point guassian quadrature
GQ.point = zeros(3,2);
GQ.weight = zeros(3,1);
GQ.point(1,1) = 0.5; GQ.point(1,2) = 0;
GQ.point(2,1) = 0;   GQ.point(2,2) = 0.5;
GQ.point(3,1) = 0.5;  GQ.point(3,2) = 0.5;
GQ.weight(1:3) = 1/6;

% set up the shape function
phi = zeros(3,3); % shape function
dphi = zeros(2,3,3); % ksi and eta derivative of shape function
for i = 1: 3 % number of GQ point
    ksi = GQ.point(i,1);
    eta = GQ.point(i,2);
    % shape function for tiangular element
    phi(1,i) = 1- ksi - eta;
    phi(2,i) = ksi;
    phi(3,i) = eta;
    % xi derivative of s
    dphi(1,1,i) = -1;
    dphi(1,2,i) = 1;
    dphi(1,3,i) = 0;
    % eta derivative of s
    dphi(2,1,i) = -1;
    dphi(2,2,i) = 0;
    dphi(2,3,i) = 1;
    
end

% assembly of stiffness and force matrices
K = zeros(numnod,numnod);
F = zeros(numnod,1);
M = zeros(numnod,numnod);

% loop over each element
for e = 1:numele
    
    % Translate shape derivative to wrt xi and eta
    e_coord = zeros(3,2); % node coordinate in global system
    for j = 1: 3 % 3 node triangular shape
        index = node(j,e);
        e_coord(j,1) = x(index);
        e_coord(j,2) = y(index);
    end
    
    % for each GQ point
    for k = 1:3
        Jcob(:,:) = dphi(:,:,k)*e_coord(:,:);
        % derivative  of shape functions wrt x and y at each gq point
        dphi_xy (:,:,k)= inv(Jcob(:,:))*dphi(:,:,k);
        % determinant of the jacobian wrt at xi and eta
        detJcob(k) = det(Jcob);
    end
    
    % calculate element stiffness matrix and mass matrix
    [ke,me,fe] = elementstiff_diff(node,e,phi,dphi_xy,detJcob,k_B,D_T,GQ,x,y);
    
    % assembly the global matrix
    for i = 1:3
        index = node(i,e);
        for j = 1:3
            index2 = node(j,e);
            K(index,index2)= K(index,index2)+ke(i,j);
            M(index,index2)= M(index,index2)+me(i,j);
        end
    end
    
    for i = 1: 3
        index = node(i,e);
        F(index)=  F(index)+ fe(i);
    end
    
    
end
A = zeros(numnod,numnod);
A = M-dt*K;

ifixu = zeros(1,numnod);
dbar = zeros(1,numnod);
for i = 1: ndivl+1
    ifixu(i) = 1; 
    ifixu(i+(ndivw*(ndivl+1)))= 2;
end 

fext_o = zeros(numnod,1);
fext = zeros(numnod,1);
% apply dirichlet boundary conditions at the edge 
for i = 1:numnod % number of node in a row 
    % set the boundary condition of the nodes in row 1 
    if(ifixu ==1)
        for j = 1: numnod 
            fext_o(j) = fext_o(j) - A(j,i)*dbar(n);
        end 
    A(i,:) = zeros(1,numnod);
    A(:,i) = zeros(numnod,1);
    A(i,i) = c0;
    end 
     if(ifixu ==2)
    A(i,:) = zeros(1,numnod);
    A(:,i) = zeros(numnod,1);
    end 
end
invA = inv(A);

% loop over all time steps
% apply backward Euler 
Tspan = 7200; 
Nsteps = Tspan/dt;
time = 0; 

d = zeros(1,numnod);
v = zeros(1,numnod);
d_tilde = zeros(1,numnod);
% initialize the time step 
for i = 1:numnod
    d(i) = c0;
end 
% time loop over 2 hr 
for i = 1: Nsteps
    time = time + dt; 
    fext = fext_o; 
    d_tilde = d; 
    % find the known f term on the right side of equation
    fext = fext - K*d_tilde';
    for j = 1: numnod
        if(ifixu(j)==1)
            fext(j) = dbar(j);
        end 
    end 
    % calculate v n+1 
    v = (invA*fext)';
    % find d n+1 
    d = d_tilde + dt*v;
    % store d n+1 
    soln_d (i,:) = d;
   
end

% dmap1 = reshape(soln_d(2,:,:),8,6);
% dmap1 = dmap1';
% surf(dmap1)
% title('C(x,y,t) at t = 1s')
% xlabel('length (node)')
% ylabel('heigth (node)')
% zlabel('Concentration (g/ml)')
% 
% figure
% dmap3600 = reshape(soln_d(3600,:,:),8,6);
% dmap3600 = dmap3600';
% surf(dmap3600)
% title('C(x,y,t) at t = 1hr')
% xlabel('length (node)')
% ylabel('heigth (node)')
% zlabel('Concentration (g/ml)')
% 
% figure
% dmap7200 = reshape(soln_d(7200,:,:),8,6);
% dmap7200 = dmap7200';
% surf(dmap7200)
% title('C(x,y,t) at t = 2hr')
% xlabel('length (node)')
% ylabel('heigth (node)')
% zlabel('Concentration (g/ml)')

%% Testing Method of Manufacturing solution 
tspan = linspace(0,Tspan,Nsteps);
for i = 1: Nsteps 
    dm(i,:) = c0+(1/tspan(i)).*exp(x.*y); 
end 
dm(1,:)= c0;
fm = SourceFun(x,y,dm,D_T,k_B,tspan,c0); 

% assembly of stiffness and force matrices
K_m = zeros(numnod,numnod);
M_m = zeros(numnod,numnod);

% loop over each element
for e = 1:numele
    
    % Translate shape derivative to wrt xi and eta
    e_coord = zeros(3,2); % node coordinate in global system
    for j = 1: 3 % 3 node triangular shape
        index = node(j,e);
        e_coord(j,1) = x(index);
        e_coord(j,2) = y(index);
    end
    
    % for each GQ point
    for k = 1:3
        Jcob(:,:) = dphi(:,:,k)*e_coord(:,:);
        % derivative  of shape functions wrt x and y at each gq point
        dphi_xy (:,:,k)= inv(Jcob(:,:))*dphi(:,:,k);
        % determinant of the jacobian wrt at xi and eta
        detJcob(k) = det(Jcob);
    end
    
    % calculate element stiffness matrix and mass matrix
    [ke,me,fe] = elementstiff_diffManuSol(node,e,phi,dphi_xy,detJcob,D_T,GQ,x,y,fm,k_B);
    
    % assembly the global matrix
    for i = 1:3
        index = node(i,e);
        for j = 1:3
            index2 = node(j,e);
            K_m(index,index2)= K_m(index,index2)+ke(i,j);
            M_m(index,index2)= M_m(index,index2)+me(i,j);
        end
    end
    
    
end
A_m = zeros(numnod,numnod);
A_m = M_m-dt*K_m;
fext_om = zeros(numnod,1);
fext_m = zeros(numnod,1);
% apply dirichlet boundary conditions at the edge 
for i = 1:numnod % number of node in a row 
    % set the boundary condition of the nodes in row 1 
    if(ifixu ==1)
        for j = 1: numnod 
            fext_om(j) = fm(2,j);
        end 
    A_m(i,:) = zeros(1,numnod);
    A_m(:,i) = zeros(numnod,1);
    A_m(i,i) = c0;
    end 
     if(ifixu ==2)
    A_m(i,:) = zeros(1,numnod);
    A_m(:,i) = zeros(numnod,1);
    end 
end
invA_m = inv(A_m);

d = zeros(1,numnod);
v = zeros(1,numnod);
d_tildem = zeros(1,numnod);
% initialize the time step 
for i = 1:numnod
    d(i) = c0;
end 
% time loop over 2 hr 
for i = 1: Nsteps
    time = time + dt;
    fext_m = fext_om ; 
    d_tildem = d; 
    % find the known f term on the right side of equation
    fext_m = fext_m - K_m*d_tildem';
    for j = 1: numnod
        if(ifixu(j)==1)
            fext_m(j) = dbar(j);
        end 
    end 
    % calculate v n+1 
    v = (invA_m*fext_m)';
    % find d n+1 
    d = d_tildem + dt*v;
    % store d n+1 
    soln_dm (i,:) = d;
end

% figure 
% plot(dm(3600,:))
% hold on 
% plot(soln_dm(3600,:))
% title('Compare exact manufactured solution with fem solution ')
% xlabel('node')
% ylabel('Concentration')
% legend('Manufactured solution','fem solution')