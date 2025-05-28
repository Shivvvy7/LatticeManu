function Main(file_path)
    clc;
    disp('=== MATLAB SCRIPT STARTED ===');

    clearvars -except file_path

    %% Parameters (match old script)
    nx = 50;
    ny = 50;
    nz = 50;

    outputDir = 'outputs';
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    %% Read STL and voxelize (same as old)
    [OUTPUTgrid1] = VOXELISE(nx, ny, nz, file_path, 'xyz');
    OUTPUTgrid = OUTPUTgrid1;

    %% Check if voxel grid is non-empty
    if nnz(OUTPUTgrid) == 0
        error('Voxelized grid is empty. Check STL file or resolution.');
    end

    %% Element counts
    nelx = size(OUTPUTgrid, 2);
    nelz = size(OUTPUTgrid, 3);
    nely = size(OUTPUTgrid, 1);
    Slab_thickness = nely;

    %% Thermal FEA initialization
    KE_th = lk_H8_ther(1);
    
    nele_th = nelx * Slab_thickness * nelz;
    ndof_th = (nelx + 1) * (Slab_thickness + 1) * (nelz + 1);
    nodegrd_th = reshape(1:(Slab_thickness + 1)*(nelx + 1), Slab_thickness + 1, nelx + 1);
    nodeids_th = reshape(nodegrd_th(1:end - 1, 1:end - 1), Slab_thickness * nelx, 1);
    nodeidz_th = 0:(Slab_thickness + 1)*(nelx + 1):(nelz - 1)*(Slab_thickness + 1)*(nelx + 1);
    nodeids_th = repmat(nodeids_th, size(nodeidz_th)) + repmat(nodeidz_th, size(nodeids_th));
    edofVec_th = nodeids_th(:) + 1;

    offsets = [0, Slab_thickness + 1, Slab_thickness, -1, ...
               (Slab_thickness + 1)*(nelx + 1), ...
               (Slab_thickness + 1)*(nelx + 2), ...
               Slab_thickness + (Slab_thickness + 1)*(nelx + 1), ...
               -1 + (Slab_thickness + 1)*(nelx + 1)];

    edofMat_th = repmat(edofVec_th, 1, 8) + repmat(offsets, nele_th, 1);
    iKt_th = reshape(kron(edofMat_th, ones(8,1))', [], 1);
    jKt_th = reshape(kron(edofMat_th, ones(1,8))', [], 1);

    %% Loading and boundary conditions
    top_vec = 1:Slab_thickness + 1:(nelx * (Slab_thickness + 1)) + 1;
    bottom_vec = Slab_thickness + 1:Slab_thickness + 1:(nelx + 1)*(Slab_thickness + 1);
    no_ele_xy = (Slab_thickness + 1)*(nelx + 1);
    q_in = top_vec;
    fixeddofs = bottom_vec;

    for itr = 1:nelz
        q_in = [q_in, top_vec + itr * no_ele_xy];
        fixeddofs = [fixeddofs, bottom_vec + itr * no_ele_xy];
    end

    %% Thermal simulation
    k_min = 1e-10; q_min = 1e-10;
    con_penal = 9; con_flux = 9;
    conductivity = 42; per_ele = 50;
use_xPhys = OUTPUTgrid;
    Temp_field = slab_analysis3D_vec(OUTPUTgrid, iKt_th, jKt_th, KE_th, ...
        conductivity, k_min, Slab_thickness, nelx, nelz, ...
        nele_th, q_in, fixeddofs, con_penal, per_ele, q_min, con_flux);

    %% Density VTK output
    [xx, yy, zz] = meshgrid(1:nx, 1:ny, 1:nz);
    new_map = ones(ny, nx, nz);
    new_map(OUTPUTgrid < 0.5) = 0;

    vtkwrite(fullfile(outputDir, 'density.vtk'), 'structured_grid', ...
             xx, yy, zz, 'scalars', 'map', new_map);

    %% Temperature VTK output
    buffer = (1/8) * ones(2,2,2);
    OUTPUTgrid_nodal = convn(OUTPUTgrid, buffer);  % remove 'same' to match old

    [xxn, yyn, zzn] = meshgrid(1:nx+1, 1:ny+1, 1:nz+1);
    Temp_field(OUTPUTgrid_nodal < 0.5) = 0;  % match threshold from old script

    vtkwrite(fullfile(outputDir, 'Temp_out.vtk'), 'structured_grid', ...
             xxn, yyn, zzn, 'scalars', 'temp', Temp_field);

    %% Temperature stats — exact same as old code
    T_top = Temp_field(:,:,end);

   T_avg = mean(T_top(T_top > 0));
T_max = max(T_top(T_top > 0));
   T_avg_3d = sum(Temp_field(:)) / numel(Temp_field);

    disp(['Average Temperature of Topmost Layer: ', num2str(T_avg)]);
    disp(['Maximum Temperature of Topmost Layer: ', num2str(T_max)]);
    disp(['Average Temperature of the Whole Lattice: ', num2str(T_avg_3d)]);
    disp('=== MATLAB SCRIPT ENDED ===');
end
