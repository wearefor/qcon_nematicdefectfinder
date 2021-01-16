function phi=func_defectorient(centroids,ring,params,nx,ny,defectstr)

% Defect Orientation Finder called by func_defectfind()
% ------------------------------------------------------------
% Michael M. Norton, Physics @ Brandeis Univeristy, 2017-2021
% ------------------------------------------------------------
% Approach: Draw a circle around a defect core and look for the region on
% that loop with a director field that most closely aligns with the radial
% unit vector of that circle.
% For +1/2 defects we report the propulsion direction for "active extensile nematics,"
% which is this vector + pi
% ------------------------------------------------------------
% [map,map_p,map_m,centroids_p,centroids_m,phi_p,phi_m]=func_defectfind(params,nx,ny)
% ------------------------------------------------------------
% inputs:
%      1. nx,ny : director
%      2. defectstr : 'p' for +1/2, 'm' for -1/2 defects 
%      3. ring : recycle the path used for line integral, but something
%      else could be used here
%
% outputs: 
%      1. phi : angle in radians measured from the x-axis
% ------------------------------------------------------------   

[x1,y1]=find(ring==1);
x1=x1-(ceil(params.N_window/2)+2);
y1=y1-(ceil(params.N_window/2)+2);
    
N_def=size(centroids,1);
phi=zeros(N_def,1);

[H,W]=size(nx);

for i=1:N_def
    
    x0=centroids(i,1);
    y0=centroids(i,2);
    
    x=x0+x1;
    y=y0+y1;
    
    % handle defects near boundary
    x(x>W)=W;
    x(x<1)=1;
    y(y>H)=H;
    y(y<1)=1;
    
   theta=atan2(y-y0,x-x0);
   theta=theta+2*pi*(theta<0);
   [theta,theta_i]=sort(theta);
    
    r1=ceil([x,y]);
    r1=r1(theta_i,:); %sorts r1 from 0 to 2*pi
    
    ny_local=ny(sub2ind(size(ny),r1(:,2)',r1(:,1)'));
    nx_local=nx(sub2ind(size(nx),r1(:,2)',r1(:,1)'));

    dotprod=abs(nx_local.*cos(theta')+ny_local.*sin(theta'));
    [~,mindot]=max(dotprod);
    
    switch defectstr
        case 'p'
            phi(i)=theta(mindot)+pi;
        case 'm'
            phi(i)=theta(mindot);    
    end
    
end