function [map,map_p,map_m,centroids_p,centroids_m,phi_p,phi_m] = func_defectfind(params,nx,ny)

% Defect Location and Orientation Finder
% ------------------------------------------------------------
% Michael M. Norton, Physics @ Brandeis Univeristy, 2017-2021
% ------------------------------------------------------------
% Approach: Use Matlab's conv2() to calculate line
% integrals on the nematic orientation field. Defect free regions will
% have a value ~0, regions containing a net topological charge will have a
% positive or negative value to around ~+/-0.5 (this threshold can be set).
% Defects are then identified using regionprops().
% ------------------------------------------------------------
% [map,map_p,map_m,centroids_p,centroids_m,phi_p,phi_m]=func_defectfind(params,nx,ny)
% ------------------------------------------------------------
% inputs:
%      1. params.N_window : radius of line integral
%      2. params.defectthresh : threshold for , something around 0.1-0.2
%
% outputs: 
%      1. map : charge map
%      2,3. map_p(m) : map showing only positive(negative) regions used to
%         identify +(-) 1/2 defects
%      4,5. centroids_p(m) : location of postive(negative) defects
%      6,7. phi_p(m) : orientation of positive(negative) defects in radians, measured from x-axis
% ------------------------------------------------------------    

[filter_x,filter_y,ring] = func_gencircle(ceil(params.N_window/2),params.gencircdebug);

Qxx=nx.^2-1/2;
Qxy=nx.*ny;

[Qxx_x,Qxx_y]=gradient(Qxx);
[Qxy_x,Qxy_y]=gradient(Qxy);

dphidx=2*(-2*Qxy.*Qxx_x+(1+2*Qxx).*Qxy_x)./(1+4*Qxx+4*Qxy.^2+4*Qxx.^2);
dphidy=2*(-2*Qxy.*Qxx_y+(1+2*Qxx).*Qxy_y)./(1+4*Qxx+4*Qxy.^2+4*Qxx.^2);

%apply filters to gradients to calculate line integral:
% Int (dphidy . dy + dphidx . dx))
% map is a matrix of the local winding number
map=(conv2(dphidy,filter_y,'same')-conv2(dphidx,filter_x,'same'))/(2*pi);

%map(abs(map)>1)=0; %threshold to eliminate spurious values
%map_p=map>params.defectthresh; %threshold for positive defects
%map_m=map<-params.defectthresh; %threshold for negative defects

dphidx=2*(-2*Qxy.*Qxx_x+(1+2*Qxx).*Qxy_x)./(1+4*Qxx+4*Qxy.^2+4*Qxx.^2);
dphidy=2*(-2*Qxy.*Qxx_y+(1+2*Qxx).*Qxy_y)./(1+4*Qxx+4*Qxy.^2+4*Qxx.^2);

num_dphidx=2*(-2*Qxy.*Qxx_x+(1+2*Qxx).*Qxy_x);
num_dphidy=2*(-2*Qxy.*Qxx_y+(1+2*Qxx).*Qxy_y);
denom=(1+4*Qxx+4*Qxy.^2+4*Qxx.^2);

eps_mine=1e0;

dphidx_filt=dphidx;
dphidy_filt=dphidy;

dphidx_filt(abs(denom)<eps_mine & abs(num_dphidx)<eps_mine) =0;
dphidy_filt(abs(denom)<eps_mine & abs(num_dphidy)<eps_mine) =0;
map_filt=(conv2(dphidy_filt,filter_y,'same')-conv2(dphidx_filt,filter_x,'same'))/(2*pi);

map=map_filt;

map_p=map>params.defectthresh; %threshold for positive defects
map_m=map<-params.defectthresh; %threshold for negative defects


%% some cleanup of the maps, not always needed

se = strel('disk',1);
%imerode(originalBW,se);

map_p=imerode(map_p,se);
map_m=imerode(map_m,se);
map_p=imdilate(map_p,se);
map_m=imdilate(map_m,se);
map_p=imclose(map_p,se);
map_m=imclose(map_m,se);


%% search maps for defects, this could be improved a bit.. 
% MinorAxisLength and Area are not dimensionless so should be adjusted

s_p = regionprops(map_p,'centroid','MinorAxisLength','Area','Extent','Eccentricity');
s_m = regionprops(map_m,'centroid','MinorAxisLength','Area','Extent','Eccentricity');

centroids_p = cat(1, s_p.Centroid);
centroids_m = cat(1, s_m.Centroid);

ext_p = cat(1, s_p.Extent);
ext_m = cat(1, s_m.Extent);

% regions containing isolated defects will generally be circular,
% artifacts can therefore be identified by aspect ratio
%minoraxislength_p = cat(1, s_p.MinorAxisLength);
%minoraxislength_m = cat(1, s_m.MinorAxisLength);

%ecc_p = cat(1, s_p.Eccentricity);
%ecc_m = cat(1, s_m.Eccentricity);

areas_p = cat(1, s_p.Area);
areas_m = cat(1, s_m.Area);


%select centroids based on some criteria
%centroids_p=centroids_p(minoraxislength_p>2 & areas_p>10 & ext_p>0.35 & ecc_p < 0.9,:);
%centroids_m=centroids_m(minoraxislength_m>2 & areas_m>10 & ext_m>0.4 & ecc_m < 0.9,:);

centroids_p=centroids_p(areas_p>10 & ext_p>0.35,:);
centroids_m=centroids_m(areas_m>10 & ext_m>0.4,:);


%% find defect orientations
%[~,~,ring] = func_gencircle(ceil(params.N_window/2),0);

[phi_p]=func_defectorient(centroids_p,ring,params,nx,ny,'p');
[phi_m]=func_defectorient(centroids_m,ring,params,nx,ny,'m');

%quiver(x0,y0,cos(theta_sort(mindot)+pi),sin(theta_sort(mindot)+pi),10,'m','LineWidth',3)

   
end

