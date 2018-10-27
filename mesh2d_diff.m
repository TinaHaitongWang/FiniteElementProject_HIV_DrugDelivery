% CEE 530 Final Project 
% Author: Haitong Wang 
% 2D FEM Diffusion problem (time-dependent)

% mesh geenration program for 2d diffusion 
function [x,y,node,numele, numnod,ndivl,ndivw] = mesh2d_diff(length,height)

% number of element in each direction 
ndivl = 7; % x 
ndivw = 5; % y

numele = ndivw*ndivl*2; 
numnod = (ndivl+1)*(ndivw+1);

% set up nodal coordinates
for i = 1: (ndivl+1)
    for j = 1: (ndivw+1)
        x((ndivl+1)*(j-1)+i) = (length/ndivl)*(i-1); 
        y((ndivl+1)*(j-1)+i) = (height/ndivw)*(j-1);
    end 
end 
   

% set up connectivity array: for 3 node triangular element 
elem = 0; 
nodet = zeros(numele,3);
for i = 1: numele
        elem = elem +1;
        if elem ==1 
            nodet (elem, 1) = elem ;
            nodet (elem, 2) = elem +1;
            nodet (elem, 3) = ndivl +1 +1;
        elseif mod(elem,2) == 0 % element is even
            if mod(elem, 7) == 0 % right bound element 
                nodet (elem, 1) = nodet (elem-1, 2);
                nodet (elem, 2) = nodet (elem-1, 3) +1;
                nodet (elem, 3) =  nodet (elem-1, 3);
            else % element is even
                nodet (elem, 1) = nodet (elem-1, 1) +1;
                nodet (elem, 2) = nodet (elem-1, 3) +1;
                nodet (elem, 3) =  nodet (elem-1, 3);
            end 
        elseif mod(elem,2) == 1 % element is odd
             if mod(elem, 7) == 1 % left bound element 
                nodet (elem, 1) = nodet (elem-1, 1)+1;
                nodet (elem, 2) = nodet (elem-1, 1)+2;
                nodet (elem, 3) =  nodet (elem-1, 2)+1;
            else % element is odd
                nodet (elem, 1) = nodet (elem-1, 1);
                nodet (elem, 2) = nodet (elem, 1) +1;
                nodet (elem, 3) =  nodet (elem-1, 2);
            end 
        
      
        end 
        
end 
        
node = nodet';

end 