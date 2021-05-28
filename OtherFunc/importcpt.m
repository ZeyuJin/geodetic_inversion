function cmap=importcpt(fname)
  % cmap=IMPORTCPT(fname) computes a MATLAB-style colormap from the GMT-
  % compatible .cpt RGB-file fname created by makecpt. Values specified
  % by 'B' and 'F' (background and foreground) will be honored, but 'N'
  % (NaN) specifications are ignored.
  %
  % Example:
  %   cmap=importcpt('~/Desktop/custom_jet.cpt');
  %   colormap(cmap);
  %
  % Eric Lindsey, Nov. 2011
  
  %read the file. This may fail if you have a highly customized cpt file.
  fid=fopen(fname);
  cmapcell=textscan(fid,'%f %d %d %d %f %d %d %d','commentstyle','#');
  Bcell=textscan(fid,'B %d %d %d','commentstyle','#');
  Fcell=textscan(fid,'F %d %d %d','commentstyle','#');
  fclose(fid);
  
  %convert 0-255 integer colormap to double (range 0-1)
  cmap=double([cmapcell{2} cmapcell{3} cmapcell{4}; ...
        cmapcell{6}(end) cmapcell{7}(end) cmapcell{8}(end)])/255;
      
  %add background and foreground colors, if specified and unique
  B=[Bcell{1} Bcell{2} Bcell{3}];
  if(~isempty(B))
    dB=double(B)/255;
    if(dB ~= cmap(1,:))
      cmap=[dB; cmap];
    end
  end
  
  F=[Fcell{1} Fcell{2} Fcell{3}];
  if(~isempty(F))
    dF=double(F)/255;
    if(dF ~= cmap(end,:))
      cmap=[cmap; dF];
    end
  end    
  
end