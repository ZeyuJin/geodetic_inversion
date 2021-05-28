function plot_errorbar(xin,yin,sig_y,varargin);
% plot the errobar 
%
%
%%%%%% Available options %%%%%%%%%%%%%
% LineWidth
% MarkerSize
% MakkerEdgeColor
% MarkerFaceColor
% LineColor
% BarWidth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Kang Wang on 11/30/2016


% default values
LineWidth=1;
MarkerSize=5;
MarkerFaceColor='r';
MarkerEdgeColor='k';
LineColor='k';
BarWidth=0.01;

if ~isempty(varargin);
    for c=1:floor(length(varargin)/2);
      try
          switch varargin{c*2-1};
              case 'LineWidth'
                  LineWidth=real(varargin{c*2});
              case 'MarkerSize'
                  MarkerSize=real(varargin{c*2});
              case 'MarkerFaceColor'
                  MarkerFaceColor=varargin{c*2};
              case 'MarkerEdgeColor'
                  MarkerEdgeColor=varargin{c*2};
              case 'LineColor'
                  LineColor=varargin{c*2};
              case 'BarWidth'
                  BarWidth=real(varargin{c*2});
          end
      catch
          error(['Unrecognized Keywords: ',varargin{c*2-1},'\n']);
    end
    end

end    


X=xin;
Y=yin;
E=sig_y;

xlength=BarWidth*(max(X)-min(X));
%xlength=0.02;
plot([X, X]',[Y-E, Y+E]',[LineColor,'-'],'LineWidth',LineWidth);
hold on
plot(X,Y,'o','MarkerSize',MarkerSize,'MarkerFaceColor',MarkerFaceColor,'MarkerEdgeColor',MarkerEdgeColor);

for k=1:length(X);
  x = [X(k) - xlength, X(k) + xlength];
  y_h = [Y(k) + E(k), Y(k) + E(k)];
  line(x, y_h,'LineWidth',LineWidth,'Color',LineColor);
  y_b = [Y(k) - E(k), Y(k) - E(k)];
  line(x, y_b,'LineWidth',LineWidth,'Color',LineColor);
end
hold off
