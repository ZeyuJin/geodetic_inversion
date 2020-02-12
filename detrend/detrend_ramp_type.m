function [ramp_type,add_col] = detrend_ramp_type(ramp_index)
   switch ramp_index
       case 4
           ramp_type = 'bi_ramp';
           add_col = 4;
       case 5
           ramp_type = 'qu_ramp_5';
           add_col = 5;
       case 7
           ramp_type = 'qu_ramp_7';
           add_col = 7;
       case 0
           ramp_type = 'noramp';
           add_col = 0;
       otherwise
           disp('There is something wrong with the index number!');
           return;
   end
end