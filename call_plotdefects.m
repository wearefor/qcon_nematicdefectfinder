% ------------------------------------------------------------
% qconv_nematicdefectfinder
% Michael M. Norton, Physics @ Brandeis Univeristy, 2017-2021
% ------------------------------------------------------------
% This script shows an example of how to use the defect finder functions.
% The example uses polscope data, which provides a direct measure of the
% nematic orientation field {nx,ny}, but you can use other methods to
% identify the orientation field. 

% info about example data:
% experimental data acquired by Linnea Lemma @ U.C. Santa Barbara
% https://pubs.rsc.org/en/content/articlelanding/2020/sm/d0sm01316a/

% (i)   orientation :: local nematic orientation represented as an angle
%       black = 0, white = 180 degrees --> in a nematic theta=theta+180 <--
%       images acquired using polscope imaging, see Pol-Acquisition-Manual-2.0.pdf for details
% (ii)  fluorescence :: microtubules labeled with alexa647 dye

%
% Note the discontinuities in the orientation field - these can be a pain
% for both plotting the nematic field and any post-processing needs that
% require differentiation
%
% (a) plotting tricks:
% to avoid "parts" in quiver plots, plot both the director n and -n:
    %quiver(nx,ny); quiver(-nx,-ny)
%    
% (b) introduce the quadratic nematic order tensor
            % Q=s|nx.^2-1/2    nx*ny  | = |Qxx Qxy|
            %    |  nx*ny    ny.^2-1/2|   |Qyx Qyy|
            %
            % some properties:
            %
            % nx^2+ny^2=1
            %
            % Q(-n)=Q(n)
            %
            % Qxy=Qyx, Qxx=-Qyy
            %
            % s=2*sqrt(Qxx^2+Qxy^2)
%
% ------------------------------------------------------------
% ------------------------------------------------------------

close all
clear all
clc


%% load data

dir='~/sampleframe';
%dir='change'
orientation=imread(['sampleframe/polscope.tif']);
%retardance=imread([dir 'img_000000000_1_Retardance - Computed Image_000.tif']);
fluorescence=imread(['sampleframe/fluorescence.tif']);

%% some pre-processing of nematic director
%[h,w]=size(orientation);
params.gauss=4; %amount of blurring to director

% define director
theta=double(orientation)*pi/(18000);
nx=cos(theta);
ny=-sin(theta);

% define Q tensor and nematic order S
Qxx=nx.^2-1/2;
Qxy=nx.*ny;

% note:
% since {nx,ny} are discontinuous they should never be smooth as this will create artifacts,
% smooth {Qxx,Qxy} representation 
Qxx_smooth=imgaussfilt(Qxx,params.gauss); 
Qxy_smooth=imgaussfilt(Qxy,params.gauss);
S=2*sqrt((Qxx_smooth.^2+Qxy_smooth.^2)); %nematic order parameter

%normalize Q
Qxx_norm=Qxx_smooth./S; 
Qxy_norm=Qxy_smooth./S; 

%redefine n
nx_smooth=sqrt(Qxx_norm+1/2);       
ny_smooth=sqrt(1-nx_smooth.^2).*sign(Qxy_norm);

%overwrite
nx=nx_smooth;
ny=ny_smooth;


% find defects, see help on func_defectfind_20190630() for details
params.N_window=10;         %size of loop integral
params.defectthresh=0.15;   %charge threshold, plot the "map" to get an idea of what a good value is for other images
params.gencircdebug=0;      %1 to show line integral filter
tic
[map,map_p,map_m,centroids_p,centroids_m,phi_p,phi_m] = func_defectfind(params,nx,ny);
toc

%% plot

fig5=figure('NumberTitle', 'on', 'Name', 'Nematic Order S=2*sqrt(Qxx^2+Qxy^2), Director, Defects');
%imshow(imadjust(fluorescence))
%imagesc(S)
imagesc(imadjust(fluorescence));
colormap(fig5,gray)
colorbar


hold on
colorbar

% to help with plotting, throw out regions where {nx,ny} flips to {-nx,-ny}
theta2=atan2(ny,nx);
[theta_x,theta_y]=gradient(theta2);
theta_grad_mag=theta_x.^2+theta_y.^2;
nx_plot=nx;
ny_plot=ny;
nx_plot(theta_grad_mag>0.05)=NaN;
ny_plot(theta_grad_mag>0.05)=NaN;


h_slice=streamslice(nx_plot,ny_plot,12,'noarrows');
set(h_slice,'Color','k','LineWidth',1);

% plot defects
defect_arrow_scale=0.1;
plot(centroids_p(:,1),centroids_p(:,2),'o','MarkerFaceColor','m','MarkerEdgeColor','k')
plot(centroids_m(:,1),centroids_m(:,2),'o','MarkerFaceColor','c','MarkerEdgeColor','k')

%plot defect arrows
quiver(centroids_p(:,1),centroids_p(:,2),cos(phi_p),sin(phi_p),defect_arrow_scale,'m','LineWidth',3)
quiver(centroids_m(:,1),centroids_m(:,2),cos(phi_m),sin(phi_m),defect_arrow_scale,'.c','LineWidth',3)
quiver(centroids_m(:,1),centroids_m(:,2),cos(phi_m+2*pi/3),sin(phi_m+2*pi/3),defect_arrow_scale,'.c','LineWidth',3)
quiver(centroids_m(:,1),centroids_m(:,2),cos(phi_m+4*pi/3),sin(phi_m+4*pi/3),defect_arrow_scale,'.c','LineWidth',3)

axis square


