function [ph_new, in, out] = mask_phase(x,y,ph,mask_area,varargin)
% this function is only used to mask the phase at near field
% x,y should be the matrix form (using meshgrid function)
% mask_area must be the closed polygon
   xv = mask_area(:,1)';
   yv = mask_area(:,2)';
   [in,on] = inpolygon(x,y,xv,yv);
   out = ~inpolygon(x,y,xv,yv);
   
   % the third option 'all' will skip all steps behind
   if nargin > 4
      for arg = 1:nargin-4
          if strcmpi(varargin{arg},'allout')
             ph_new = ph;
             ph_new(in) = NaN;
             ph_new(on) = NaN;
             return;
          elseif strcmpi(varargin{arg},'allin')
             ph_new = ph;
             ph_new(out) = NaN;
             return;
          end
      end
   end

   % mask_area is the two_column array of vertex
   out = ~inpolygon(x,y,xv,yv);
   ph1 = ph; ph1(in) = 0;
   ph2 = ph; ph2(out) = 0;
   
   % check whether it is positive or negative area
   ph3 = ph(in); tmp = nanmean(ph3);
   if tmp > 0
       ph2(ph2 < 0) = NaN;
   elseif tmp < 0
       ph2(ph2 > 0) = NaN;
   else
       warning('Unexpected masking, exit!');
       return;
   end

   ph_new = ph1 + ph2;
   ph_new(on) = NaN; % prevent some points on the edge of polygon
end
