 function [ status ] = is_leapyear( year )
if mod(year, 400) == 0
status = 1;
elseif mod(year, 4) == 0 && mod(year, 100) ~= 0
status = 1;
else
status = 0;
end
end