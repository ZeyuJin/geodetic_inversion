function out = errorbarxy(x,y,varargin)
%ERRORBARXY Customizable error bar plot in X and Y direction
%
%   This function allows the user to plot the graph of x against y, along with
%   both x and y errorbars. 
% 
%   With 4 numeric arguments (x,y,dx,dy), error bar are assumed to be of 
%   same magnitude in both direction.
%
%   One can specify lower and upper error bar with 6 numeric arguments
%   (x,y,dx_high,dy_high,dx_low,dy_low).
%
%   x,y,dx,dy,... must be vectors of the same length
%
%   It is possible to customize the line properties of the error bars by
%   adding pair of 'field/value' fields (such as 'LineWidth',2) that can be
%   understood by line. See LineProperties for more information. Entering
%   simple string like 'ko' will be ignored.
%
%   OUTPUT
%   h = errorbarxy(...) returns h, a vector of 2 graphic handles, one for
%   the line between points, one for the error bars.
% 
%   EXAMPLES
%   X = 10 * rand(7,1);
%   Y = 10 * rand(7,1);
%   dx = rand(7,1);
%   dy = rand(7,1);
%   errorbarxy(X,Y,dx,dy,'Color','k','LineStyle','none','Marker','o',...
%   'MarkerFaceColor','w','LineWidth',1,'MarkerSize',11);
%   figure
%   X = 10 * rand(7,1);
%   Y = 10 * rand(7,1);
%   dx = rand(7,1);
%   dy = rand(7,1);
%   dx2 = rand(7,1);
%   dy2 = rand(7,1);
%   errorbarxy(X,Y,dx,dy,dx2,dy2,'Color','B','LineStyle','--','Marker','s',...
%   'MarkerFaceColor','w','LineWidth',2,'MarkerSize',11);
%
%   ACKNOWLEDGEMENT
%   This is a rewrite of the m-file errorbarxy of James Rooney, to add
%   customizable line properties.

% Version info:
%
% - Nov 2008: from Keith Brady feedback, prevent oddities when entering
% negative error bar length
%
% ------------------ INFO ------------------
%   Authors: Jean-Yves Tinevez 
%   Work address: Max-Plank Institute for Cell Biology and Genetics, 
%   Dresden,  Germany.
%   Email: tinevez AT mpi-cbg DOT de
%   November 2007 - November 2008;
%   Permission is given to distribute and modify this file as long as this
%   notice remains in it. Permission is also given to write to the author
%   for any suggestion, comment, modification or usage.
% ------------------ BEGIN CODE ------------------

    %% Parse argins and set options
    
    nargs = numel(varargin);
    for argi = 1 : nargs        
        if ~( isnumeric(varargin{argi}) )
            break
        end
        errbaropt{argi} = varargin{argi};        
    end    
    displayopt = varargin(argi:end);
    
    options.Color = 'k';
    
    % Get display options
    if ~isempty(displayopt)        
        if isstruct(displayopt{1})
            options = displayopt{1};
        elseif numel(displayopt) > 1
            options = varargin2struct(displayopt);
        end        
        erroroptions = options;        
    end
    
    % For errorbar lines
    erroroptions.LineStyle = '-';
    erroroptions.Marker = 'none';
        
    % Lateral size of error bars
    xw = (max(x)-min(x))/100;
    yw = (max(y)-min(y))/100;
    
    % Get numeric input (error limits)
    n = numel(errbaropt);    
    if n == 2
        % only 2 cells, so this is the same for lower and upper bar
        ux = errbaropt{1};
        lx = ux;
        uy = errbaropt{2};
        ly = uy;        
    elseif n == 4
        % 4 cells, the user specified both upper and lower limit
        ux = errbaropt{1};
        lx = errbaropt{3};
        uy = errbaropt{2};
        ly = errbaropt{4};        
    else
        errid = 'MATLAB:errorbarxy:BadArgumentNumber';
        errmsg = ['Must have 2, 4 or 6 numeric arguments, got ' ,num2str(n+2),'.'];
        error(errid,errmsg);        
    end
    
    % To prevent oddities with negative error bars
    ux = abs(ux);
    lx = abs(lx);
    uy = abs(uy);
    ly = abs(ly);
    
    
    %% Actual plotting
    
    holdstate = get(gca,'NextPlot');
    X = [];
    Y = [];
    for t=1:length(x)
        
        % x errorbars
        X = [ X     nan x(t)-lx(t) x(t)+ux(t)    nan    x(t)-lx(t) x(t)-lx(t) nan     x(t)+ux(t) x(t)+ux(t)         ];
        Y = [ Y     nan y(t) y(t)                nan    y(t)-yw y(t)+yw       nan     y(t)-yw y(t)+yw               ];
        
        % y errorbars
        X = [ X     nan x(t) x(t)                nan    x(t)-xw x(t)+xw       nan     x(t)-xw x(t)+xw               ];
        Y = [ Y     nan y(t)-ly(t) y(t)+uy(t)    nan    y(t)-ly(t) y(t)-ly(t) nan     y(t)+uy(t) y(t)+uy(t)         ];
        
    end
    
    hp = plot(x,y,options);
    hold on
    he = plot(X,Y,erroroptions);
    
    uistack(he,'top');
    uistack(hp,'top');
    out = [hp he];
    
    % Return to initial hold state 
    set(gca,'NextPlot',holdstate)
    
    %% Subfunctions
    
    function out = varargin2struct(in)
        
        if ~iscell(in)
            errid = 'MATLAB:struct2varargin:BadInputType';
            errmsg = ['Input argument must be a cell, got a ' ,class(in),'.'];
            error(errid,errmsg);
        end
        
        n = length(in);
        
        if mod(n,2) ~= 0
            errid = 'MATLAB:struct2varargin:BadInputType';
            errmsg = ['Input argument must have an even number of elements, got ' ,num2str(n),'.'];
            error(errid,errmsg);
        end
        
        out = struct;
        
        for i = 1 : n/2
            name = in{2*i-1};
            value = in{2*i};
            out.(name) = value;
        end
    end

end
