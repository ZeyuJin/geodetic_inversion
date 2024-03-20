function read_plot_fault_segment(filename, color, option, varargin)

    % option is used to specify what registration to be used
    % lon/lat or utm x/y

   %addpath('/Users/zej011/work/Kang_tutorial/codes_utilities/matlab/igppsar');
   
   data = load(filename);
   lon = data(:,1);
   lat = data(:,2);
   
   % default parameters
   lon_eq = 99;
   lat_eq = 34;
   ref_lon = 99;
   if ~isempty(varargin)
       for CC = 1:floor(length(varargin)/2)
           try
               switch lower(varargin{CC*2-1})
                   case 'lonc'
                       lon_eq = varargin{CC*2};
                   case 'latc'
                       lat_eq = varargin{CC*2};
                   case 'ref_lon'
                       ref_lon = varargin{CC*2};
               end
           catch
               error('Unrecognized Keyword');
           end
       end
   end
   [xo,yo] = ll2xy(lon_eq,lat_eq,ref_lon);
   
   
   % different fault segments are separated by 370/0
   tmp_index = find(abs(lon - 370) < eps);
   fault_num = length(tmp_index);
   seg_index = zeros(1,(fault_num+1)); % add 0 index at the beginning for slice use
   seg_index(2:end) = tmp_index;
   clear tmp_index
   

   for i=1:fault_num
       j1 = seg_index(i) + 1;    % start index of fault segment
       j2 = seg_index(i+1) - 1;  % end index of fault segment

       % shrink the number of data points
%        p = plot(lon(j1:j2),lat(j1:j2),'b','linewidth',3);
%        tmp_fx = round(fx * ratio);
%        tmp_fy = round(fy * ratio);
%        [~,idx] = unique([tmp_fx,tmp_fy],'rows');
%        count = count + length(idx);
%        disp(length(idx));

       if strcmp(option, 'll')
           px = lon(j1:j2);
           py = lat(j1:j2);
       else
           [fx,fy] = ll2xy(lon(j1:j2),lat(j1:j2),ref_lon);
           px = (fx - xo) / 1000;
           py = (fy - yo) / 1000;
       end
       
       plot(px,py,'color',color,'linewidth',2);
       
   end
end
