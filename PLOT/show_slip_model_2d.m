function show_slip_model_2d(slip_model,varargin)
% Plotting the 2-d profile of the slip mode
% Xiaoyu Zou 04/09/2025
format long
%set(0,'defaultAxesFontName', 'AvantGarde')
%set(0,'defaultAxesFontSize', 15)
d2r=pi/180;

% data=load(slip_model);
data = slip_model;
% iflt=data(:,1);
% indx=data(:,2);
% indx_layer=data(:,3);
xe=data(:,4);
yn=data(:,5);
zr=data(:,6);
lp=data(:,7);
wp=data(:,8);
p_strk=data(:,9);
p_dip=data(:,10)*d2r;
tp=data(:,11);
% slip1=-data(:,12)/100;
% slip2=-data(:,13)/100; %May need to adjust it based on which side of the block you are looking at.
% min_s1=min(slip1);
% max_s1=max(slip1);
% min_s2=min(slip2);
% max_s2=max(slip2);

% nflt=length(unique(iflt));


% W=sum(unique(wp));
% Ny=length(unique(indx_layer));
% N=length(xe);

cmax = 500;
seis_matfile = [];
% axis_range = [-30 20 0 50 -25 0];
fault_file = [];
lon_eq = -117.5; 
lat_eq = 35.5;
ref_lon = lon_eq;
% title = 'cm';

if ~isempty(varargin)
   for CC = 1:floor(length(varargin)/2)
       try
          switch lower(varargin{CC*2-1})
              case 'misfit_range'
                  cmax = varargin{CC*2};
              case 'seismic'
                  seis_matfile = varargin{CC*2};
              case 'axis_range'
                  axis_range = varargin{CC*2};
              case 'fault_file'
                  fault_file = varargin{CC*2};
              case 'ref_lon'
                  ref_lon = varargin{CC*2};
              case 'lonc'
                  lon_eq = varargin{CC*2}; 
              case 'latc'
                  lat_eq = varargin{CC*2};
%               case 'title'
%                   title = varargin{CC*2};
          end
       catch
          error('Unrecognized Keyword');
       end
   end
end




for k=1:max(slip_model(:,1))
    this_seg=slip_model(find(slip_model(:,1)==k),:);
    ddl=unique(this_seg(:,8));%down dip length of this segment

    figure()
    hold on
    SLIP=[];
    XV=[];
    YV=[];
    ZV=[];
    XO=[];
    YO=[];
    ZO=[];
    for n=1:max(this_seg(:,3))
        this_layer=this_seg(find(this_seg(:,3)==n),:);
        if n==1
            xstart=0;ystart=0;
        else
            xstart=0;ystart=ystart-ddl(n-1);
        end
        
        N=size(this_layer,1);
        lp=this_layer(:,7);
        wp=this_layer(:,8);
        slip1=this_layer(:,12)/100;
        slip2=this_layer(:,13)/100;
        for i=1:N
            if i==1
                x1f=xstart;y1f=ystart;
            else
                x1f=x2f;y1f=y2f;
            end
            

            x2f=x1f+lp(i);
            y2f=y1f;


            x3f=x2f;
            y3f=y2f-wp(i);


            x4f=x1f;
            y4f=y3f;




            X=[x1f,x2f,x3f,x4f]/1000;
            Y=[y1f,y2f,y3f,y4f]/1000;


            xo=mean(X);
            yo=mean(Y);

            C=sqrt(slip1(i)^2+slip2(i)^2);
            patch(X,Y,C,'FaceAlpha',1);
            hold on
            xuf=slip1(i);
            yvf=slip2(i);





            XV=[XV;xuf];
            YV=[YV;yvf];


            XO=[XO;xo];
            YO=[YO;yo];


            slip0=sqrt(xuf^2+yvf^2);
            SLIP=[SLIP;slip0];

        end
        slipmax=max(SLIP);
        quiver(XO,YO,XV/slipmax,YV/slipmax,1.5,'LineWidth',0.5,'Color',[0,0,0]);
    end


    if ~isempty(seis_matfile)
        dseis = load(seis_matfile);
        slon = dseis(:,1);
        slat = dseis(:,2);
        sdepth = dseis(:,3) * -1;

        [xo,yo] = ll2xy(lon_eq,lat_eq,ref_lon);
        [xx,yy] = ll2xy(slon,slat,ref_lon);
        xs = (xx - xo) ./ 1000;
        ys = (yy - yo) ./ 1000;

        xv = [-3 -5 20 20 -2];
        yv = [33 51 51 30 32];
        out = ~inpolygon(xs,ys,xv,yv);
        xin = xs(out);
        yin = ys(out);
        din = sdepth(out);

        S = 5;
        scatter3(xin,yin,din,S,'black','filled');
    end

    % if ~isempty(fault_file)
    %     for ii = 1:LS
    %         slon = [lonf(ii) lonf(ii+LS)];
    %         slat = [latf(ii) latf(ii+LS)];
    %         [xx,yy] = ll2xy(slon,slat,ref_lon);
    %         xs = (xx - xo) ./ 1000;
    %         ys = (yy - yo) ./ 1000;
    %         line(xs,ys,'color','black','linewidth',1.5);
    %     end
    % end

    xlabel('Along-Strike (km)');
    ylabel('Along-Dip (km)');

    axis equal
    % zlim([-25 0]);
    % axis([-60 -10 285 305 -25 0]);
    % axis(axis_range);
    colormap jet
    % hc = colorbar('southoutside');
    % hc = colorbar('southoutside');
    hc = colorbar;
    title(hc,'slip (m)');
    if cmax < 5
        caxis([-cmax cmax]);
    else
        caxis([0 cmax/100]);
    end
    % colorbar off
    % oldcmap = colormap;
    % colormap( flipud(oldcmap) );

    grid on
    set(gca,'Fontsize',20,'fontweight','bold');
    set(gcf,'PaperPositionMode','auto');
end

end
