function ori = mesh_norm(mesh)
% Compute normals from the coregistred brain mesh using the patchnormals
% function
% If patchnormals fails for any triangles, the Fieldtrip function
% "normals" is used
%
%-CREx180530

mesh.vertices = mesh.pos;
mesh.faces = mesh.tri;
ori = patchnormals(mesh);
isn = isnan(ori(:, 1));
if any(isn)
    warning('Patchnormals function failed to find normals for %s vertex - Using fieldtrip normals',...
                num2str(sum(isn)));
    ftnrm = normals(mesh.pos, mesh.tri);
    ori(isn==1, :) = ftnrm(isn==1, :);
end