function H_okada=calc_smooth_okada(slip_model_in);
%generate the smooth matrix given the slip model
%by Kang Wang on 03/16/2015
%Last Updated by Kang Wang on 03/16/2015
%Last Updated by Kang Wang on 08/20/2015
format long
d2r=pi/180.0;

%slip_model=load_PRM(PRMFILE,'SLIP_MODEL_IN');
%data=load(slip_model);
data=slip_model_in;
indx=data(:,2);
indx_layer=data(:,3);
xe=data(:,4);
yn=data(:,5);
zr=data(:,6);
lp=data(:,7);
wp=data(:,8);
p_strk=data(:,9);
p_dip=data(:,10)*d2r;

Np=length(indx);

xfo=zeros(Np,1);
yfo=zeros(Np,1);
zfo=zeros(Np,1);
xo=zeros(Np,1);
yo=zeros(Np,1);
zo=zeros(Np,1);

sp=lp.*wp;   %area of patches
%sp=ones(Np,1);

%Ny=max(indx_layer); %number of layers of the patches
Ny=length(unique(indx_layer)); %number of layers

strike_a=mean(p_strk);
theta_x=(90.0-strike_a)*d2r;
p_theta=zeros(size(p_strk));

for i=1:Np
    p_theta(i)=(90.0-p_strk(i))*d2r;
    [x1f,y1f]=xy2XY(xe(i),yn(i),p_theta(i));
    z1f=zr(i);
        
    x2f=x1f+lp(i);
    y2f=y1f;
    z2f=z1f;
    
    x3f=x2f;
    y3f=y2f-wp(i)*cos(p_dip(i));
    z3f=z2f-wp(i)*sin(p_dip(i));
    
    x4f=x1f;
    y4f=y3f;
    z4f=z3f;
    xfo(i)=mean([x1f,x2f,x3f,x4f]);
    yfo(i)=mean([y1f,y2f,y3f,y4f]);
    zfo(i)=mean([z1f,z2f,z3f,z4f]);
end

[xeo,yno]=xy2XY(xfo,yfo,-p_theta);
[xo,yo]=xy2XY(xeo,yno,theta_x);
zo = zfo;   % added by Zeyu Jin

H1=zeros(Np,Np); %smooth matrix for strike slip
H2=zeros(Np,Np); 

for i=1:Ny;
   if (i>1&i<Ny);
   indx_ith=find(indx_layer==i);
   indx_top=find(indx_layer==(i-1));
   indx_bottom=find(indx_layer==(i+1));
  
   Np_ith=length(indx_ith);
   Np_top=length(indx_top);
   Np_bottom=length(indx_bottom);
   
   for j=1:Np_ith;
    %   indx_patch=length(find(indx_layer<i))+j;      
       indx_patch=indx_ith(j);
       
       x1=xo(indx_ith(j));
       y1=yo(indx_ith(j));
       z1=zo(indx_ith(j));
     
       dis_top=zeros(Np_top,1);
       for k=1:Np_top;
           x2=xo(indx_top(k));
           y2=yo(indx_top(k));
           z2=zo(indx_top(k));
           dis_top(k)=p2p_dis([x1,y1,z1],[x2,y2,z2]);
       end
       [dis1,n1]=sort(dis_top);
       indxt=indx_top(n1(1));
  
       dis_bottom=zeros(Np_bottom,1);       
       for l=1:Np_bottom;
           x3=xo(indx_bottom(l));
           y3=yo(indx_bottom(l));
           z3=zo(indx_bottom(l));
           dis_bottom(l)=p2p_dis([x1,y1,z1],[x3,y3,z3]);
       end
       [dis2,n2]=sort(dis_bottom);
       indxb=indx_bottom(n2(1));
       
       if (j==1);
           H1(indx_patch,indx_ith(j+1))=-1;
%            H1(indx_patch,indxt)=-1*sp(indxt)/sp(indx_patch+1);
%            H1(indx_patch,indxb)=-1*sp(indxb)/sp(indx_patch+1);
           H1(indx_patch,indxt)=-1;
           H1(indx_patch,indxb)=-1;
           H1(indx_patch,indx_patch)=-(H1(indx_patch,indx_patch+1)+H1(indx_patch,indxt)+H1(indx_patch,indxb));
       elseif(j==Np_ith);
           H1(indx_patch,indx_ith(j-1))=-1;
%            H1(indx_patch,indxt)=-1*sp(indxt)/sp(indx_patch-1);
%            H1(indx_patch,indxb)=-1*sp(indxb)/sp(indx_patch-1);
           H1(indx_patch,indxt)=-1;
           H1(indx_patch,indxb)=-1;
           H1(indx_patch,indx_patch)=-(H1(indx_patch,indx_patch-1)+H1(indx_patch,indxt)+H1(indx_patch,indxb));
       else
           H1(indx_patch,indx_ith(j-1))=-1;
           H1(indx_patch,indx_ith(j+1))=-1;
%            H1(indx_patch,indxt)=-1*sp(indxt)/sp(indx_patch-1);
%            H1(indx_patch,indxb)=-1*sp(indxb)/sp(indx_patch-1);
           H1(indx_patch,indxt)=-1;
           H1(indx_patch,indxb)=-1;
           H1(indx_patch,indx_patch)=-(H1(indx_patch,indx_patch-1)+H1(indx_patch,indx_patch+1)+H1(indx_patch,indxt)+H1(indx_patch,indxb));
       end
       
   end
    
  elseif (i==1); 
          indx_ith=find(indx_layer==1);
          indx_bottom=find(indx_layer==2);
         
          Np_ith=length(indx_ith);
          Np_bottom=length(indx_bottom);
          
          for j=1:Np_ith;
       	      indx_patch=indx_ith(j);
              
              x1=xo(indx_ith(j));
              y1=yo(indx_ith(j));
              z1=zo(indx_ith(j));
              dis_bottom=zeros(Np_bottom,1);       
              for l=1:Np_bottom;
                 x3=xo(indx_bottom(l));
                 y3=yo(indx_bottom(l));
                 z3=zo(indx_bottom(l));
                  dis_bottom(l)=p2p_dis([x1,y1,z1],[x3,y3,z3]);
              end
              [dis2,n2]=sort(dis_bottom);
              indxb=indx_bottom(n2(1));
              
              if (j==1);
                  H1(indx_patch,indx_ith(j+1))=-1;
%                   H1(indx_patch,indxb)=-1*sp(indxb)/sp(indx_patch+1);
                  H1(indx_patch,indxb)=-1;
                  H1(indx_patch,indx_patch)=-(H1(indx_patch,indx_patch+1)+H1(indx_patch,indxb));
              elseif (j==Np_ith);
                  H1(indx_patch,indx_ith(j-1))=-1;
%                   H1(indx_patch,indxb)=-1*sp(indxb)/sp(indx_patch-1);
                  H1(indx_patch,indxb)=-1;
                  H1(indx_patch,indx_patch)=-(H1(indx_patch,indx_patch-1)+H1(indx_patch,indxb));
              else
                  H1(indx_patch,indx_ith(j-1))=-1;
                  H1(indx_patch,indx_ith(j+1))=-1;
%                   H1(indx_patch,indxb)=-1*sp(indxb)/sp(indx_patch+1);
                  H1(indx_patch,indxb)=-1;
                  H1(indx_patch,indx_patch)=-(H1(indx_patch,indx_patch-1)+H1(indx_patch,indx_patch+1)+H1(indx_patch,indxb));
              end   
          end
   else
          indx_ith=find(indx_layer==Ny);
          indx_top=find(indx_layer==Ny-1);
          
          Np_ith=length(indx_ith);
          Np_top=length(indx_top);
          
          for j=1:Np_ith;
              %indx_patch=length(find(indx_layer<Ny))+j;
                     indx_patch=indx_ith(j);

              x1=xo(indx_ith(j));
              y1=yo(indx_ith(j));
              z1=zo(indx_ith(j));
              
              dis_top=zeros(Np_top,1);
              for k=1:Np_top;
                  x2=xo(indx_top(k));
                  y2=yo(indx_top(k));
                  z2=zo(indx_top(k));
                  dis_top(k)=p2p_dis([x1,y1,z1],[x2,y2,z2]);
              end
              [dis1,n1]=sort(dis_top);
              indxt=indx_top(n1(1));
              if (j==1);
                H1(indx_patch,indx_ith(j+1))=-1; 
%                 H1(indx_patch,indxt)=-1*sp(indxt)/sp(indx_patch+1);
                H1(indx_patch,indxt)=-1;
                H1(indx_patch,indx_patch)=-(H1(indx_patch,indx_patch+1)+H1(indx_patch,indxt));
              elseif(j==Np_ith);
                H1(indx_patch,indx_ith(j-1))=-1; 
%                 H1(indx_patch,indxt)=-1*sp(indxt)/sp(indx_patch-1);
                H1(indx_patch,indxt)=-1;
                H1(indx_patch,indx_patch)=-(H1(indx_patch,indx_patch-1)+H1(indx_patch,indxt));        
              else
                H1(indx_patch,indx_ith(j-1))=-1;
                H1(indx_patch,indx_ith(j+1))=-1;
%                 H1(indx_patch,indxt)=-1*sp(indxt)/sp(indx_patch+1);
                H1(indx_patch,indxt)=-1;
                H1(indx_patch,indx_patch)=-(H1(indx_patch,indx_patch-1)+H1(indx_patch,indx_patch+1)+H1(indx_patch,indxt));
              end
          end
   end

end
   
FDip = 5; % dip-slip smoothing is heavier than strike-slip 
H2 = FDip .* H1;
H_okada=[H1,zeros(size(H1));zeros(size(H2)),H2];
end
