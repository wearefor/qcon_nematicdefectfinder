function [filter_x,filter_y,ring] = func_gencircle(rad,debug)

% Integral Path, called by func_defectfind()
% ------------------------------------------------------------
% Michael M. Norton, Physics @ Brandeis Univeristy, 2017-2019
% ------------------------------------------------------------
% Makes a pixelated ring for calculating line integrals.
% This can probably be simplified.. but the script only needs to be run
% once to make the filter.
% ------------------------------------------------------------
% [filter_x,filter_y,ring]=func_gencircle(rad,debug)
% inputs:
%      1. rad : radius of ring
%      2. debug : 1 plots the ring with the tangent vector, 0 does nothing
%
% outputs: 
%      1. filter_x(y) : coordinates along ring
%      2. ring : ring "mask" / logical
% ------------------------------------------------------------


%generates a ring of pixels
disk_outer=fspecial('disk',rad+1);
N=length(disk_outer);

disk_outer(disk_outer==disk_outer(floor(rad),floor(rad)))=1;
disk_outer(disk_outer<disk_outer(floor(rad),floor(rad)))=0;


disk_inner=fspecial('disk',rad);
disk_inner(disk_inner==disk_inner(floor(rad),floor(rad)))=1;
disk_inner(disk_inner<disk_inner(floor(rad),floor(rad)))=0;


pad_size=(length(disk_outer)-length(disk_inner))/2;
disk_inner_big=padarray(disk_inner,[pad_size pad_size],0,'both');

disk_outer(disk_inner_big==1)=0;
ring=disk_outer;
%imagesc(disk_outer)
%pause
[x_grid,y_grid]=meshgrid((1:N)-N/2,(1:N)-N/2);

r_grid=sqrt(x_grid.^2+y_grid.^2);

x_grid(disk_outer==0)=0;
y_grid(disk_outer==0)=0;

if debug==1
    size(disk_outer)
    imagesc(disk_outer); hold on;
    quiver(-y_grid./r_grid,x_grid./r_grid,0.5,'m')
    pause
    close all
else
end

filter_x=-y_grid./r_grid;
filter_y=-x_grid./r_grid;

end

