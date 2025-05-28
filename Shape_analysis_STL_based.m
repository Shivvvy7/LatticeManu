clc
clear all
%% reading the STL
nx=20;
ny=20;
nz=20;
[OUTPUTgrid1] = VOXELISE(nx,ny,nz,'TPMS_Gyroid.stl','xyz');
OUTPUTgrid = OUTPUTgrid1;

%---------------------------------------------
% visualize the voxalized result (optional)
% Show the voxelised result:me
figure;
subplot(1,3,1);
imagesc(squeeze(sum(OUTPUTgrid,1)));
colormap(gray(256));
xlabel('Z-direction');
ylabel('X-direction');
axis equal tight

subplot(1,3,2);
imagesc(squeeze(sum(OUTPUTgrid,2)));
colormap(gray(256));
xlabel('Z-direction');
ylabel('Y-direction');
axis equal tight

subplot(1,3,3);
imagesc(squeeze(sum(OUTPUTgrid,3)));
colormap(gray(256));
xlabel('Y-direction');
ylabel('X-direction');
axis equal tight

%% Getting number of elements
nelx=size(OUTPUTgrid,2);
nelz=size(OUTPUTgrid,3);
nely=size(OUTPUTgrid,1);
Slab_thickness = nely;
%% Thermal FEA initialization
KE_th = lk_H8_ther(1);
nele_th=nelx*Slab_thickness*nelz;
ndof_th = (nelx+1)*(Slab_thickness+1)*(nelz+1);
nodegrd_th = reshape(1:(Slab_thickness+1)*(nelx+1),Slab_thickness+1,nelx+1);
nodeids_th = reshape(nodegrd_th(1:end-1,1:end-1),Slab_thickness*nelx,1);
nodeidz_th = 0:(Slab_thickness+1)*(nelx+1):(nelz-1)*(Slab_thickness+1)*(nelx+1);
nodeids_th = repmat(nodeids_th,size(nodeidz_th))+repmat(nodeidz_th,size(nodeids_th));
edofVec_th = nodeids_th(:)+1;
edofMat_th = repmat(edofVec_th,1,8)+...
    repmat([0 Slab_thickness+1 Slab_thickness -1 (Slab_thickness+1)*(nelx+1) (Slab_thickness+1)*(nelx+2) Slab_thickness+((Slab_thickness+1)*(nelx+1)) -1+((Slab_thickness+1)*(nelx+1))],nele_th,1); 
iKt_th = reshape(kron(edofMat_th,ones(8,1))',8*8*nele_th,1);
jKt_th = reshape(kron(edofMat_th,ones(1,8))',8*8*nele_th,1);
%% Analysis initialization
nele=nelx*Slab_thickness*nelz;
per_ele=50;
k_min=1e-10;
q_min=1e-10;
con_penal=9;
con_flux=9;
conductivity=42;
%% THARMAL LOADING FOR EACH SLAB
top_vec=1:Slab_thickness+1:(nelx*(Slab_thickness+1))+1;
bottom_vec=Slab_thickness+1:Slab_thickness+1:(nelx+1)*(Slab_thickness+1);
no_ele_xy=(Slab_thickness+1)*(nelx+1);
q_in=top_vec;
fixeddofs=bottom_vec;
for itr=1:1:nelz
    store1=top_vec+(itr*no_ele_xy);
    store2=bottom_vec+(itr*no_ele_xy);
    q_in=[q_in store1];
    fixeddofs=[fixeddofs store2];
end
%% Thermal ANALYSIS
use_xPhys = OUTPUTgrid;
Temp_field =slab_analysis3D_vec(use_xPhys,iKt_th,...
       jKt_th,KE_th,conductivity,k_min,Slab_thickness,nelx,nelz,...
       nele,q_in,fixeddofs,con_penal,per_ele,q_min,con_flux);
%% writing the vtk file to be visualized in paraview
x=1:nx;
y=1:ny;
z=1:nz;
[xx,yy,zz]=meshgrid(x,y,z);


new_map=ones(ny,nx,nz);
new_map(OUTPUTgrid<0.5)=0;
filename='density.vtk';
title='map';
vtkwrite(filename,'structured_grid',xx,yy,zz,'scalars',title,new_map)

x=1:nx+1;
y=1:ny+1;
z=1:nz+1;
[xx,yy,zz]=meshgrid(x,y,z);

filename='Temp_out.vtk';
title='temp';
vtkwrite(filename,'structured_grid',xx,yy,zz,'scalars',title,Temp_field)


