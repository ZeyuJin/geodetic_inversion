function [xo,yo]=ll2xy(xi,yi,lon_c)
%
% [xo,yo]=ll2xy(xi,yi,lon_c);
%
% Last updated by Kang Wang on 11/09/2015
yold=yi;
if (yi<0)
yi=-yi;
end


%r_dtor=1.74532925199d-2;   %  = pi/180
r_dtor=pi/180;
% i_ft=0;
        i_ft=0;
        a_griddes=['C','D','E','F','G','H','J',...
            'K','L','M','N','P','Q','R','S','T','U',...
            'V','W','X'];
        r_a=6378206.4d0;
        r_e2=0.00676865799761d0;
        r_k0=0.9996d0;    %scale at center 
        r_lat0=0.0d0;
        r_fe=5e5;
        r_fn(1)=0;
        r_fn(2)=1e7;  % 
        
        r_ep2 = r_e2/(1.d0 - r_e2);
        r_e4 = r_e2^2;
        r_e6 = r_e2^3;
        r_dtor = pi/180;
        
           r_v=zeros(length(xi(:)),2);
           xi=xi*r_dtor;
           yi=yi*r_dtor;

           i_zone = fix(mod(xi+3.d0*pi,2.d0*pi)/(r_dtor*6.d0)) + 1;
           i_zone = max(min(i_zone,60),1);
         % r_lon0 = -pi + 6.d0*r_dtor*(i_zone-1) + 3.d0*r_dtor % central meridian
           %r_lon0 = lon_c;
           r_lon0=lon_c*r_dtor;
           r_n = r_a./sqrt(1.d0 - r_e2*sin(yi).^2);
           r_t = tan(yi).^2;
           r_t2 = r_t.^2;
           r_c = r_ep2*cos(yi).^2;
           r_ba = (xi - r_lon0).*cos(yi);
           r_a2 = r_ba.^2;
           r_a3 = r_ba.*r_a2;
           r_a4 = r_ba.*r_a3;
           r_a5 = r_ba.*r_a4;
           r_a6 = r_ba.*r_a5;
           r_m = r_a.*((1.0d0-r_e2/4 - 3.0d0*r_e4/64.0d0-5.d0*r_e6/256.d0)* ...
                 yi - (3.d0*r_e2/8.d0 + 3.d0*r_e4/32.d0 +45.d0*r_e6/ ...
                 1024.d0)*sin(2.d0*yi) +  (15.d0*r_e4/256.d0 + ...
                 45.d0*r_e6/ 1024.d0)*sin(4.d0*yi) - (35.d0*r_e6/3072.d0)* ...
                 sin(6.d0*yi));
           r_m0 = r_a.*((1.d0-r_e2/4 - 3.d0*r_e4/64.d0 -5.d0*r_e6/256.d0)* ...
                  r_lat0 - (3.d0*r_e2/8.d0 + 3.d0*r_e4/32.d0 + ...
                  45.d0*r_e6/ 1024.d0)*sin(2.d0*r_lat0) +  (15.d0*r_e4/ ...
                  256.d0 +45.d0*r_e6/1024.d0)*sin(4.d0*r_lat0) - ...
                  (35.d0* r_e6/3072.d0)*sin(6.d0*r_lat0));

           r_v(:,1) = r_k0*r_n.*(r_ba+(1.d0-r_t+r_c).*r_a3/6.d0 ...
                    +(5.d0-18.d0*r_t+r_t2+72.d0*r_c-58.d0*r_ep2).*r_a5/120.d0);
           r_v(:,1) = r_v(:,1) + r_fe;

           r_v(:,2) = r_k0*(r_m - r_m0 + r_n.*tan(yi).*( r_a2/2.d0 + ...
                   (5.d0-r_t+ 9.d0*r_c+4.d0*r_c.^2).*(r_a4/24.d0) + ...
                   (61.d0- 58.d0*r_t+r_t2+600.d0*r_c-330.d0*r_ep2).* ...
                   (r_a6/720.d0) ));
           if yi >= 0
              r_v(:,2) = r_v(:,2) + r_fn(1);
           else
              r_v(:,2) = r_v(:,2) + r_fn(2);
           end

           r_k = r_k0*(1.d0+(1.d0+r_ep2.*cos(yi).^2).*(r_v(:,1)-r_fe).^2./ ...
                (2.d0*(r_k0.^2).*r_n.^2));

           i_gi = fix((yi./r_dtor+80.d0)/8.d0) + 1;
           i_gi = max(min(i_gi,20),1);
           a_grid = a_griddes(i_gi);

           xo=r_v(:,1);
         if (yold<0) 
           yo=-r_v(:,2);
         else
           yo=r_v(:,2);
         end

end