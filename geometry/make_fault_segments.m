function slip_model = make_fault_segments(fault_id,xstart,ystart,dz_start,strike,dip,L,W,N_layer,lp_top,bias_lp,bias_wp)
% generate a vertical fault model given the starting point and fault
% dimension
%
% Usage:
% slip_model=make_fault_segments(fault_id,xstart,ystart,strike,L,W,N_layer,lp_top,l_bias,w_bias);
%

d2r=pi/180;

% strike
theta=(90-strike)*d2r;
[xo_segment,yo_segment]=xy2XY(xstart,ystart,theta);

wp_factor=zeros(N_layer,1);
for k=1:N_layer
  wp_factor(k)=bias_wp^(k-1);
end
wp_top=W/sum(wp_factor);

wp_layer=zeros(N_layer,1);
for k=1:N_layer
  wp_layer(k)=wp_top*bias_wp^(k-1);
end

indx_patch=0;   %% it will be recomputed finally to combine all the fault segments
% dz_start=0;   % fault maybe in depth

slip_model=[];
for j=1:N_layer
   lp_this_layer_rough=lp_top*bias_lp^(j-1);
   N_this_layer=round(L/lp_this_layer_rough);
   if N_this_layer == 0
       lp_this_layer = L;     % in case that one segment is really small
       N_this_layer = 1;
   else
       lp_this_layer = L/N_this_layer;
   end

   indx_fault_this_layer=zeros(N_this_layer,1);
   indx_patch_this_layer=zeros(N_this_layer,1);
   indx_depth_this_layer=zeros(N_this_layer,1);

   xp_this_layer=zeros(N_this_layer,1);
   yp_this_layer=zeros(N_this_layer,1);
   zp_this_layer=zeros(N_this_layer,1);

   wpatch_this_layer=zeros(N_this_layer,1);
   lpatch_this_layer=zeros(N_this_layer,1);
   strkp_this_layer=strike*ones(N_this_layer,1);
   dip0_this_layer=dip*ones(N_this_layer,1);
   tp_this_layer=zeros(N_this_layer,1);    % assume zero topography
   slip1_this_layer=zeros(N_this_layer,1);
   slip2_this_layer=zeros(N_this_layer,1);
   
   for k=1:N_this_layer
       indx_patch=indx_patch+1;
       indx_fault_this_layer(k)=fault_id;
       indx_patch_this_layer(k)=indx_patch;
       indx_depth_this_layer(k)=j;
       xpatch_tmp=xo_segment+(k-1)*lp_this_layer;
       ypatch_tmp=yo_segment-cosd(dip)*(sum(wp_layer(1:j))-wp_layer(j));
       [xp_this_layer(k),yp_this_layer(k)]=xy2XY(xpatch_tmp,ypatch_tmp,-(theta));
       zp_this_layer(k)=(-(sum(wp_layer(1:j))-wp_layer(j)))*sind(dip)+dz_start;
       lpatch_this_layer(k)=lp_this_layer;
       wpatch_this_layer(k)=wp_layer(j);
       
       slip_model=[slip_model;indx_fault_this_layer(k),indx_patch,indx_depth_this_layer(k),xp_this_layer(k),yp_this_layer(k),zp_this_layer(k),...
                            lpatch_this_layer(k),wpatch_this_layer(k),strkp_this_layer(k),dip0_this_layer(k),tp_this_layer(k),slip1_this_layer(k),slip2_this_layer(k)];
       

   end
end