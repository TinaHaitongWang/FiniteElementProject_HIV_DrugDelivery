% CEE 530 Final Project
% Author: Haitong Wang
% 2D FEM Diffusion problem (time-dependent)

% element stiffness matrix for 3-nodes triangular element using 2-points
% gaussian quadrature

function [ke,me,fe] = elementstiff_diff(node,e,phi,dphixy,detJcob,k_B,D_T,GQ,x,y)
% loop over all gaussian point
ke = zeros(3,3);
me = zeros(3,3);
fe = zeros(3,1);

for a = 1: 3
    % find the global x and y coordinates wrt xi and eta
    xstore = 0;
    ystore = 0;
    for k = 1:3 % loop over all three nodes
        index = node(k,e);
        xstore = xstore + phi(k,a) * x(index);
        ystore = ystore + phi(k,a) * y(index);
    end
    
    % element stiffness and mass matrix
    for i = 1: 3 % loop over all three nodes
        for j = 1: 3 %
            ke(i,j) = ke(i,j) + (D_T * (dphixy(1,i,a)*dphixy(1,j,a) +...
                dphixy(2,i,a)*dphixy(2,j,a)) + k_B * phi(i,a) * phi(j,a))...
                * detJcob(a) * GQ.weight(a); 
            
            me(i,j) = me(i,j) + (phi(i,a) * phi(j,a))* detJcob(a)* GQ.weight(a);
            
            fe(i,j) = 0;
        end
    end
    
    
    
    
end


end


