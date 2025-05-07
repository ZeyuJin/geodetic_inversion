clear
% Mask out noisy values in InSAR data
addpath(genpath('/Volumes/T7/Research/PamirProject/geodetic_inversion-master'));

%% Configuration Here
path='./ALOS/data/rng/'; %path to data directory
%path = './Postseismic/DES121/';
%path='./ALOS/data/';
filepath=[path,'rng_ll_masked.grd']; %path to the original data
%maskpath=[path,'azo_ll_masked.grd'];% Path to the masked data
maskpath=[path,'rng_ll_masked.grd'];% Path to the masked data or any data for viewing only
running_mode=2;% 1: running the masking; 2: viewing the saved result
sign_mask=0;%1: do the sign masking; 0: no need to do the sign masking
load /Volumes/T7/Research/Tingri_Project/geodetic_data/Google_Earth_Data/fault_main.txt
load /Volumes/T7/Research/Tingri_Project/geodetic_data/Google_Earth_Data/fault_sub.txt

%% Running here

[X,Y,Z]=grdread2(filepath);
[x,y]=meshgrid(X,Y);
x=x(:);
y=y(:);
Z=Z(:);

if sign_mask==1 && running_mode==1
    

    figure()
    scatter(x,y,80,Z,'filled')
    hold on
    plot(fault_main(:,1),fault_main(:,2),'-ro','LineWidth',5)
    clim([-max(Z)*0.25,max(Z)*0.25])
    hold on
    title('sign-masking; masking out negative values');

    while true
        h = impoly;  % Draw a polygon interactively
        if isempty(h)
            break;
        end
        pos = getPosition(h);
        x_poly = pos(:,1);
        y_poly = pos(:,2);

        % Check which points are inside the polygon
        in = inpolygon(x, y, x_poly, y_poly);

        % Set points inside the polygon to NaN
        Z(in & Z<0)=NaN;

        %Continue drawing more polygons or break the loop
        answer = input('Do you want to draw another polygon to mask out negative values? (y/n): ', 's');
        if lower(answer) ~= 'y'
            break;
        end
    end
end


if sign_mask==1 && running_mode==1
    figure()
    scatter(x,y,80,Z,'filled')
    hold on
    plot(fault_main(:,1),fault_main(:,2),'-ro','LineWidth',5)
    clim([-max(Z)*0.25,max(Z)*0.25])
    hold on
    title('sign-masking; masking out positive values');

    while true
        h = impoly;  % Draw a polygon interactively
        if isempty(h)
            break;
        end
        pos = getPosition(h);
        x_poly = pos(:,1);
        y_poly = pos(:,2);

        % Check which points are inside the polygon
        in = inpolygon(x, y, x_poly, y_poly);

        % Set points inside the polygon to NaN
        Z(in & Z>0)=NaN;

        %Continue drawing more polygons or break the loop
        answer = input('Do you want to draw another polygon to mask out positive values? (y/n): ', 's');
        if lower(answer) ~= 'y'
            break;
        end
    end
end




if running_mode==1
    figure()
    scatter(x,y,80,Z,'filled')
    hold on
    plot(fault_main(:,1),fault_main(:,2),'-ro','LineWidth',5)
    clim([-max(Z)*0.25,max(Z)*0.25])
    hold on
    title('Draw polygons and double-click to confirm each, enter n to exit');

    while true
        h = impoly;  % Draw a polygon interactively
        if isempty(h)
            break;
        end
        pos = getPosition(h);
        x_poly = pos(:,1);
        y_poly = pos(:,2);

        % Check which points are inside the polygon
        in = inpolygon(x, y, x_poly, y_poly);

        % Set points inside the polygon to NaN
        Z(in)=NaN;

        %Continue drawing more polygons or break the loop
        answer = input('Do you want to draw another polygon? (y/n): ', 's');
        if lower(answer) ~= 'y'
            break;
        end
    end


    %% If you want to save your data:
    figure();
    scatter(x,y,80,Z,'filled')
    hold on
    plot(fault_main(:,1),fault_main(:,2),'-ro','LineWidth',5)
    clim([-max(Z)*0.25,max(Z)*0.25])
    title('Updated Scatter Plot with NaN values inside polygons')
    Z=reshape(Z,[size(Y,2),size(X,2)]);

    answer = input('Do you want to save your data? (y/n): ', 's');
    if lower(answer) == 'y'
        grdwrite2(X,Y,Z,maskpath)
    end
end



%% Second running mode:
if running_mode==2
    [X,Y,Z]=grdread2(maskpath);
    [x,y]=meshgrid(X,Y);
    x=x(:);
    y=y(:);
    Z=Z(:);

    figure()
    title('masked data')
    scatter(x,y,80,Z,'filled')
    hold on
    plot(fault_main(:,1),fault_main(:,2),'-ro','LineWidth',5)
    plot(fault_sub(:,1),fault_sub(:,2),'-ro','LineWidth',5)
    clim([-max(Z)*0.1,max(Z)*0.1])
    set(gca,'FontSize',18)
    xlabel('Lon')
    ylabel('Lat')
end

