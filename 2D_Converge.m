%% Finite Difference Problem
% 
% This code is developed to study a 2d steady state thermal problem of the form:
% 
%              Adiabatic
% 
%            __________________________
%           |                          |
%           |                          | 
%           |                          |
% Adiabatic |                          | Tinf hinf    
%           |                          | 
%           |                          |
%           |                          | 
%           |                          |
%           |                          | 
%           ---------------------------- Tb
% 
%                        
% 
% Using the finite difference method. A mesh convergance study can also be done 
% in order to understand how the number of nodes affects the accuracy of the solution.

clear

% 1-Inputs

%% Material Constants
k0 = 10;


%%Problem constants
Lx = 20; %Length x
Ly = 20; %Length y

%%NOTE: This code assumes that dx=dy. You have to make sure that this assumption is accurate
%when defining Lx, Ly, Nx, Ny.
Nx =20; %Nodes in x
Ny =20; %Nodes in y



% Dirichlet condition boundary condition on bottom
Tb= 400;
%Convective boundary condition (Robin) on right
Tinf = 1700;
hinf=1000;
%The left and top edges of the plate have adiabatic B.C.

%Spatial resolution of nodes
dx=Lx/(Nx-1);
dy=Ly/(Ny-1);

%% Define Boundary Conditions


%Assign a mapping matrix which asssigns a tag to every node.
k = 1;
node_id=zeros(Nx,Ny);
for iy = 1:Nx
    for ix=1:Ny
      node_id(ix,iy)=k;
      k=k+1;
    end
end


%A is the coefficient matrix which has dimensions (NxN) where N is the total number of nodes.
%The ith row of matrix A represents the energy balance equation for node i. Column j represents the 
%contribution of the jth node to every node. So A(i,j) is the contribution of node j to the energy balance
%at node i.
A = zeros((Nx)*(Ny)); %Create an N by N matrix and initialize with zeros

%C is the constant vector which has dimensions (N). The ith value in this vector, are the constants
%for the energy equation of the ith node.
C = zeros((Nx)*(Ny),1); %Create an N length vector and initialize with zeros

%To index this matrix and vector, you have to use the node_id as the index. Example: to access the
% node at position (2,4) you have to use A(node_id(2,4),node_id(2,4))

%Steady state solution, construct coefficient matrix assigning the
%correct forms to each node based on the node location 
for iy = 1:Ny 
    for ix=1:Nx

        if (ix == Nx) && (iy== Ny)
            %%Top right corner has 2 conductive and 2 convective terms. 
            %This node is already completed as an example.
            A(node_id(ix,iy),node_id(ix-1,iy)) = 1; %Node to the left
            A(node_id(ix,iy),node_id(ix,iy-1)) = 1; %Node to the south
            A(node_id(ix,iy),node_id(ix,iy)) = -(2 + hinf*dx/k0); %Convection surfaces
            C(node_id(ix,iy)) = -(hinf*dx/k0) * Tinf; %Constant due to convection


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fill in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %First assign the Dirichlet boundary conditions
        elseif (iy==1) 
            A(node_id(ix,iy),node_id(ix,iy)) = 1; %Convection surfaces
            C(node_id(ix,iy)) = Tb; %Constant due to convection
            

            %All nodes on the bottom edge have temp Tb, we don't have to 
            %handle the energy equations for these nodes (including the bottom 
            %left and bottom right corners).

        elseif(ix==1) &&(iy==Ny)
            
            %%Top left corner has adiabatic B.C.s
            A(node_id(ix,iy),node_id(ix+1,iy)) = 1;%Node to the right
            A(node_id(ix,iy),node_id(ix,iy-1)) = 1; %Node to the south
            A(node_id(ix,iy),node_id(ix,iy)) = -2; %Convection surfaces
            C(node_id(ix,iy)) = 0;

        elseif(ix==1)
            %%All nodes on the left edge, adiabatic boundary condition
            A(node_id(ix,iy),node_id(ix+1,iy)) = 2; %Node to the left
            A(node_id(ix,iy),node_id(ix,iy-1)) = 1; %Node to the south
            A(node_id(ix,iy),node_id(ix,iy+1)) = 1;
            A(node_id(ix,iy),node_id(ix,iy)) = -4; %Center
            C(node_id(ix,iy)) = 0; %Constant due to convection


        elseif(ix == Nx)
            %%All nodes on the right edge are convection on plane edge
            A(node_id(ix,iy),node_id(ix-1,iy)) = 2; %Node to the left
            A(node_id(ix,iy),node_id(ix,iy+1)) = 1; %Node to the north
            A(node_id(ix,iy),node_id(ix,iy-1)) = 1; %Node to the south
            A(node_id(ix,iy),node_id(ix,iy)) = -2 * (2 + hinf*dx/k0); %Convection surfaces
            C(node_id(ix,iy)) = -(2 * hinf*dx/k0) * Tinf; %Constant due to convection


            
        elseif(iy==Ny)
            %%All nodes on top edge are adiabatic on plane edge
            A(node_id(ix,iy),node_id(ix-1,iy)) = 1; %Node to the left
            A(node_id(ix,iy),node_id(ix+1,iy)) = 1;%Node to the right
            A(node_id(ix,iy),node_id(ix,iy-1)) = 2; %Node to the south
            A(node_id(ix,iy),node_id(ix,iy)) = -4; %Convection surfaces
            C(node_id(ix,iy)) = 0;

        else 
            %%All other interior nodes
            A(node_id(ix,iy),node_id(ix-1,iy)) = 1; %Node to the left
            A(node_id(ix,iy),node_id(ix+1,iy)) = 1;
            A(node_id(ix,iy),node_id(ix,iy+1)) = 1;
            A(node_id(ix,iy),node_id(ix,iy-1)) = 1; %Node to the south
            A(node_id(ix,iy),node_id(ix,iy)) = -4; %Convection surfaces
            C(node_id(ix,iy)) = 0;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

%% Now solve the matrix equation A*T=C

T = mtimes(inv(A),C);

%Now create datastructure for plotting
Tplot = zeros(Nx,Ny);
k=1;
for ix =1:Nx
    for iy=1:Ny
        Tplot(ix,iy) = T(k);
        k=k+1;
    end
end

figure

%Plot the calculated temperature field
hold on
[M,c] = contourf(0:dx:Lx,0:dy:Ly, Tplot, 20);
set(c,'edgecolor','none');
c.LineWidth= 3;

% Overlay the FDM calculation points
points=zeros(2,Nx*Ny);
i=1;
for iy=1:Ny
    for ix=1:Nx
        points(:,i) = [dx*(ix-1), dy*(iy-1)];
        i=i+1;
    end
end
hold on
plot(points(1,:),points(2,:),'r.')
hold on
rectangle('Position', [0,0,Lx,Ly])
colorbar

format long
%Interpolate and display the temperature at the center of the plate
fprintf("The center point temperature is: ")
fprintf('%0.4f',interp2(0:dx:Lx,0:dy:Ly,Tplot,Lx*0.5,Ly*0.5))