% CEE 530 Final Project
% Author: Haitong Wang
% 2D time-dependent FEM Diffusion problem

% Source term fm for manufactured solution 
function fm = SourceFun(x,y,dm,D_T,k_B,tspan,c0)
% for i = 1: length(dm(:,1))
%     fm(i,:) = (5.*x.*y.*(A-x).*(y-B))-(-10.*D_T.*tspan(i).*y.*(y-B))-(10.*D_T.*tspan(i).*x.*(A-x))...
%     + (k_B.*dm(i,:));
% end 
for i = 1: length(dm(:,1))
    fm(i,:)= (k_B.*(c0+tspan(i).*exp(x+y)))-(D_T.*tspan(i).*exp(x+y));
end 
end