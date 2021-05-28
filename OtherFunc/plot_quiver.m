function plot_quiver(x,y,u,v,varargin);
% plot the 2D vectors (and error ellipse)
% by Kang Wang on 12/01/2014
% 
% Usage: plot_quiver(x,y,u,v,'ErrorDadata',[sig_u,sig_v]);
% Available Options %%%%%%%%%%%%%%%%%
%  'Scale'
%  'ArrowColor'
%  'ArrowLineWidth'
%  'HeadLength'
%  'HeadWidth'
%  'HeadStyle' ------- {'plain', 'vback1','vback2','cback1','cback2','cback3'}
%  'ErrorData'
%  'ErrorFaceColor'
%  'ErrorEdgeColor'
%  'ErrorEdgeWidth'
%  'ErrorAlpha'

%%%%%%%%%%DEFINE THE DEFAULT VALUES
scale=1;     %scale of the vector
arrow_color='black';
arrow_line_width=1.5;
error_color='black';
error_line_color='black';
error_line_width=1.0;
head_style='vback1';
head_length=20;
head_width=10;
error_alpha=0.1;
sig_x=zeros(size(x));
sig_y=zeros(size(x));

if ~isempty(varargin);
  for c=1:floor(length(varargin)/2);
     try
        switch varargin{c*2-1};
           case 'Scale'
               scale=real(varargin{c*2});
           case 'ArrowColor'
               arrow_color=varargin{c*2};
           case 'ArrowLineWidth'
               arrow_line_width=varargin{c*2};
           case 'HeadLength'
               head_length=varargin{c*2};
           case 'HeadWidth'
               head_width=varargin{c*2};
           case 'HeadStyle'
                head_style=varargin{c*2};
           case 'ErrorData'
                error_data=varargin{c*2};
                sig_x=error_data(:,1);
                sig_y=error_data(:,2);
           case 'ErrorFaceColor'              
               error_color=varargin{c*2};
           case 'ErrorEdgeColor'
               error_line_color=varargin{c*2};
           case 'ErrorEdgeWidth'
               error_line_width=varargin{c*2};
           case 'ErrorAlpha'
               error_alpha=varargin{c*2};                 
        end
     catch
        error('Unrecognized Keyword\n');
%        fprintf('unrecognized property or value for: %s\n',varargin{c*2-1});
     end
  end
end

%x=scale*x;
%y=scale*y;
u=scale*u;
v=scale*v;
sig_x=sig_x*scale;
sig_y=sig_y*scale;
Narrow=length(x);
h0=quiver(x,y,u,v,'AutoScale','off');
set(h0,'Visible','off');
% hkid=get(h0,'children');
% X=get(hkid(1),'XData');
% Y=get(hkid(1),'YData');
X = h0.XData;
Y = h0.YData;

t = linspace(0,2*pi);
hh = zeros(size(x));
th = zeros(size(y));

dx=X(2)-X(1);
dy=Y(2)-Y(1);
L1=sqrt(dx^2+dy^2);

du=u(1);
dv=v(1);
L2=sqrt(du^2+dv^2);
scale0=L1/L2;
r=scale0;


n=0;
hold on
for i=1:3:length(X)-1;
   n=n+1;
  ah=annotation('arrow',...
      'headStyle',head_style,'HeadLength',head_length,'HeadWidth',head_width,...
      'Color',arrow_color,'LineWidth',arrow_line_width);
  set(ah,'parent',gca);
  set(ah,'position',[X(i) Y(i) X(i+1)-X(i) Y(i+1)-Y(i)]);
  xx = X(i+1) + r*sig_x(n).*cos(t).*cosd(th(n)) - r*sig_y(n).*sin(t).*sind(th(n));
  yy = Y(i+1) + r*sig_x(n).*cos(t).*sind(th(n)) + r*sig_y(n).*sin(t).*cosd(th(n));
%  hold on
  plot(xx,yy,'k-','LineWidth',0.6);
%  hp=patch(xx,yy,error_color);
%  set(hp,'FaceAlpha',error_alpha,'FaceColor',error_color,'LineWidth',error_line_width,...
%      'EdgeColor',error_line_color);
end
hold off
