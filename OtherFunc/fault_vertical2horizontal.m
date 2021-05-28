function slip_model_out=fault_vertical2horizontal(slip_model_vertical,depth_at);
% rotate a vertical fault model along its surface trace to form a horizontal fault model
% then shift downward by depth_at
%
% Usage: slip_model_out=fault_vertical2horizontal(slip_model_vertical,depth_at);
%
% by Kang Wang on 07/27/2015
% Last Updated by Kang Wang on 07/27/2015

format long
d2r=pi/180;
data0=slip_model_vertical;

indx_fault=data0(:,1);
indx_patch=data0(:,2);
indx_depth=data0(:,3);
xo_patch=data0(:,4);
yo_patch=data0(:,5);
zo_patch=data0(:,6);

lp=data0(:,7);
wp=data0(:,8);
strkp=data0(:,9);
dip0=data0(:,10);
tp=data0(:,11);
slip1=data0(:,12);
slip2=data0(:,13);

fault_id=unique(indx_fault);
Nfault=length(fault_id);

slip_model_out=[];
dip=zeros(Nfault,1);

patch_layer=unique(indx_depth);
Npatch_layer=length(patch_layer);
wp_layer=zeros(Npatch_layer,1);

for i=1:Nfault;
    indx_this_fault=find(indx_fault==fault_id(i));
    Npatch_this_fault=length(indx_this_fault);
    indx_segment=indx_fault(indx_this_fault);
    indx_patch_this_fault=indx_patch(indx_this_fault);
    indx_depth_this_fault=indx_depth(indx_this_fault);
    xo_this_fault=xo_patch(indx_this_fault);
    yo_this_fault=yo_patch(indx_this_fault);
    zo_this_fault=zo_patch(indx_this_fault);
    wp_this_fault=wp(indx_this_fault);
    lp_this_fault=lp(indx_this_fault);
    strkp_this_fault=strkp(indx_this_fault);
    dip_this_fault=dip(i)*ones(Npatch_this_fault,1);
    tp_this_fault=tp(indx_this_fault);
    slip1_this_fault=slip1(indx_this_fault);
    slip2_this_fault=slip2(indx_this_fault);
    
    theta_this_fault=(90.0-strkp_this_fault)*d2r;
    [xf_this_fault,yf_this_fault]=xy2XY(xo_this_fault,yo_this_fault,theta_this_fault);
    z_this_fault=zo_this_fault;
   
    dy_dip=zo_this_fault.*cos(dip_this_fault*d2r);
    dz_dip=zo_this_fault.*sin(dip_this_fault*d2r);
    xf_new=xf_this_fault;
    yf_new=yf_this_fault+dy_dip;
    zf_new=dz_dip;

    [xe_dip,yn_dip]=xy2XY(xf_new,yf_new,-theta_this_fault);
  data_this_fault=[indx_segment,indx_patch_this_fault,indx_depth_this_fault,xe_dip,yn_dip,zf_new,...
        lp_this_fault,wp_this_fault,strkp_this_fault,dip_this_fault,tp_this_fault,slip1_this_fault,slip2_this_fault];

  slip_model_out=[slip_model_out;data_this_fault];
end

slip_model_out(:,6)=slip_model_out(:,6)-depth_at; %shift downward

 
