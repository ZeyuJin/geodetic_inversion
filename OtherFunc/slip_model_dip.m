function slip_model_out=slip_model_dip(slip_model_in,dip_in);
%
% calculate the fault model of a given dip angle from a vertical fault model
%
% Usage:
%  slip_model_out=slip_model_dip(slip_model_in,dip_in);
%
%  slip_model_in: original (vertical) fault model
%  dip_in:  dip angle of desired fault model
%  slip_model_out: output fault model with dip angle of dip_in

% by Kang Wang on 07/26/2015
d2r=pi/180;
indx_fault=slip_model_in(:,1);
indx_patch=slip_model_in(:,2);
indx_depth=slip_model_in(:,3);

xo_patch=slip_model_in(:,4);
yo_patch=slip_model_in(:,5);
zo_patch=slip_model_in(:,6);

lp=slip_model_in(:,7);
wp=slip_model_in(:,8);
strkp=slip_model_in(:,9);
dip0=slip_model_in(:,10);
tp=slip_model_in(:,11);
slip1=slip_model_in(:,12);
slip2=slip_model_in(:,13);

fault_id=unique(indx_fault);
Nfault=length(fault_id);

Ndip=length(dip_in);
if (Nfault~=Ndip);
  error(['The # of input dipping angle must equal to # of segments!']);
end

slip_model_out=[];
dip=dip_in;

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



