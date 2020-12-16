function v = patch_smooth_filt(p,str,Nl,Nw,mode,plot_or_not)
    % For converting a okada patches solution to a spectral solution
    % p is patch paramaters that has xc,yc,zc,width,length,dip,strike
    % xs, ys and zs are the divided patch source
    % u is the corresponding slip
    % ll, ww are the corresponding length and width for the patch source
    % Nl Nw are the num of divisions of length and width for p
    
    a = load(str);
    xs = a(:,4);
    ys = a(:,5);
    zs = a(:,6);
    ll = a(:,7);
    ww = a(:,8);

    if mode == 1
        u = a(:,11);
    elseif mode == 2
        u = a(:,12);
    elseif mode == 3
        u = a(:,13);
    elseif mode == 12
        u = sqrt(a(:,11).^2+a(:,12).^2);
    elseif mode == 13
        u = sqrt(a(:,11).^2+a(:,13).^2);
    elseif mode == 23
        u = sqrt(a(:,12).^2+a(:,13).^2);
    else
        disp('ERROR: Please choose a correct mode for plotting')
        return
    end
    
    x = zeros(Nl,Nw);
    y = x;
    z = x;
    v = x;
    
    xb = [xs, xs + ll*sind(p.strike), xs + ll*sind(p.strike) + ww*cosd(p.dip)*cosd(p.strike), xs + ww*cosd(p.dip)*cosd(p.strike), xs];
    yb = [ys, ys + ll*cosd(p.strike), ys + ll*cosd(p.strike) - ww*cosd(p.dip)*sind(p.strike), ys - ww*cosd(p.dip)*sind(p.strike), ys];
    zb = [zs, zs, zs - ww*sind(p.dip), zs - ww*sind(p.dip), zs];
    
    x0 = p.x - 0.5*p.len*sind(p.strike);
    y0 = p.y - 0.5*p.len*cosd(p.strike);
    z0 = p.z;
    
    dl = p.len/Nl;
    dw = p.wid/Nw;
    
    stepxl = dl*sind(p.strike);
    stepxw = dw*cosd(p.dip)*cosd(p.strike);
    stepyl = dl*cosd(p.strike);
    stepyw = -dw*cosd(p.dip)*sind(p.strike);
    stepz = -dw*sind(p.dip);
    
    x0 = x0 + 0.5*(stepxl+stepxw);
    y0 = y0 + 0.5*(stepyl+stepyw);
    z0 = z0 + 0.5*stepz;
    
    [xx,yy]=meshgrid(1:1:Nl,1:1:Nw);
    
    x = x0 + (xx-1)*stepxl + (yy-1)*stepxw;
    y = y0 + (xx-1)*stepyl + (yy-1)*stepyw;
    z = z0 + (yy-1)*(stepz);
    %{
    x = reshape(x,numel(x),1);
    y = reshape(y,numel(y),1);
    z = reshape(z,numel(z),1);
    v = reshape(v,numel(v),1);
    %}
    
    x0 = p.x - 0.5*p.len*sind(p.strike);
    y0 = p.y - 0.5*p.len*cosd(p.strike);
    z0 = p.z;
    
    for j = 1:1:length(xs)
        bx = xb(j,:)';
        by = yb(j,:)';
        bz = zb(j,:)';
        iw1 = floor((bz(1) - z0)/stepz)+1;
        iw2 = floor((bz(3) - z0)/stepz);
        il1 = round(((bx(1) - x0) - (bz(1) - z0)/stepz*stepxw)/stepxl)+1;
        il2 = round(((bx(3) - x0) - (bz(3) - z0)/stepz*stepxw)/stepxl);
        v(iw1:iw2,il1:il2) = u(j);
    end

    h = fspecial('gaussian', 101,18);
    v = conv2(v,h,'same');
    
    %figure,surf(x,y,z,'EdgeColor','None'),hold on, plot3(xc,yc,zc,'x'),grid on,plot3(ext_x,ext_y,ext_z,'r*')
    %figure,plot3(xc,yc,zc,'o'),grid on
    %x = reshape(x,Nl*Nw,1);
    %y = reshape(y,Nl*Nw,1);
    %z = reshape(z,Nl*Nw,1);
    
    
    %v = griddata(xc,yc,zc,u,x,y,z);
    
    %x = reshape(x,Nl,Nw);
    %y = reshape(y,Nl,Nw);
    %z = reshape(z,Nl,Nw);
    %v = reshape(v,Nl,Nw);
    if plot_or_not ~= 0
        figure
        surf(x,y,z,v,'EdgeColor','None');
        colormap(jet),colorbar
        %hold on
        %plot3(xc,yc,zc,'x'),grid on,plot3(ext_x,ext_y,ext_z,'r*')
    end