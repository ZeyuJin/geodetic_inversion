function model_out=slip_model_kink(model_in,dip_range);
%
% Usage: model_out=slip_model_kink(model_in,dip_range);
%
%  modle_in: vertical slip model 
%  dip_range=[dip1,w1,w2;
%             dip2,w2,w3;
%             dip3,w3,w4];
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example:
% model_in=load('nepal_vel.slip');
% dip_range=[7,0,50e3;
%           20,50e3,120e3;
%          4,120e3,500e3];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% by Kang Wang in June, 2016

d2r=pi/180;
[m,n]=size(dip_range);

ipatch=model_in(:,2);
idepth=model_in(:,3);
xp=model_in(:,4);
yp=model_in(:,5);
zp=model_in(:,6);
lp=model_in(:,7);
wp=model_in(:,8);

indx_layer=unique(idepth);
ndepth=length(indx_layer);

npatch=length(ipatch);

wtop=abs(zp);
wbottom=abs(zp)+wp;
wcenter=abs(zp)+wp/2;

model_out=[];

dyf=0;
dw=0;
dxf=0;

indx_top=find(idepth==indx_layer(1));
ztop=zp(indx_top);
dz=mean(ztop);
for i=1:m;
   this_dip=dip_range(i,1);
   w1=dip_range(i,2);
   w2=dip_range(i,3);
   indx_this_dip=find(wcenter>=w1&wcenter<w2);
   model_vel_this_dip=model_in(indx_this_dip,:);

   wtop_this_dip=wtop(indx_this_dip);
   wbottom_this_dip=wbottom(indx_this_dip);

   xp_this_dip=model_vel_this_dip(:,4);
   yp_this_dip=model_vel_this_dip(:,5);
   zp_this_dip=model_vel_this_dip(:,6);
   lp_this_dip=model_vel_this_dip(:,7);
   wp_this_dip=model_vel_this_dip(:,8);
   strkp_this_dip=model_vel_this_dip(:,9);
   dip_this_dip=model_vel_this_dip(:,10);
   theta_this_dip=(90.0-strkp_this_dip)*d2r; 

   zp_this_dip_new=zp_this_dip+dw;
   model_vel_this_dip_new=model_vel_this_dip;
   model_vel_this_dip_new(:,6)=zp_this_dip_new;
   model_this_dip=slip_model_dip(model_vel_this_dip_new,this_dip);
   xxp_this_dip=model_this_dip(:,4);
   yyp_this_dip=model_this_dip(:,5);
   zzp_this_dip=model_this_dip(:,6);

   [xxf_this_dip,yyf_this_dip]=xy2XY(xxp_this_dip,yyp_this_dip,theta_this_dip);
   yyf_this_dip_new=yyf_this_dip+dyf;
   xxf_this_dip_new=xxf_this_dip-dxf;
   zzf_this_dip_new=zzp_this_dip;  

   [xxp_new,yyp_new]=xy2XY(xxf_this_dip_new,yyf_this_dip_new,-theta_this_dip);

   zzp_new=zzf_this_dip_new-dz;
   model_this_dip(:,4)=xxp_new;
   model_this_dip(:,5)=yyp_new;
   model_this_dip(:,6)=zzp_new;;

   model_this_dip(:,1)=i;
   model_out=[model_out;model_this_dip];
   width_this_dip=sum(unique(wp_this_dip));
   
   dxf=dxf;
   dyf=dyf-width_this_dip*cos(this_dip*d2r);
   dw=dw+width_this_dip;
   dz=dz+width_this_dip*sin(this_dip*d2r);  

end

%figure;
%show_slip_data(model_out,1);
