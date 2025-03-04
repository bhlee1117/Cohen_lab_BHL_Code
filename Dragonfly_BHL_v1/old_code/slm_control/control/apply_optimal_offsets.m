function out = apply_optimal_offsets(xyz,opt)
if ~exist('opt','var')
    opt = '';
end

z_offset = -40.5e+1;%-5.4e-1; %meter
rot_x = 60; % degree
rot_y = 0;
mat = [ % y-axis x-z plane rotation second
    cosd(rot_y) 0 sind(rot_y)
    0 1 0
    -sind(rot_y) 0 cosd(rot_y)
] * ...
[ % x-axis y-z plane rotation first
    1 0 0 
    0 cosd(rot_x)  -sind(rot_x)
    0  sind(rot_x)  cosd(rot_x)
];


switch opt
    case 'apply_mat'
        xyz_offset = mean(xyz,2);
        out = mat *(xyz-xyz_offset) + xyz_offset;
        
    case 'apply_trans'
        xyz(3,:) = xyz(3,:)+z_offset;
        out = xyz;
    case 'full'
        xyz(3,:) = xyz(3,:)+z_offset;

        xyz_offset = mean(xyz,2);
        out = mat * (xyz-xyz_offset) + xyz_offset;
    case 'mat'
        out = mat;
    case 'trans'
        out = z_offset;
end
