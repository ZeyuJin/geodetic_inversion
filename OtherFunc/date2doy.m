function output=date2doy(input);
% convert the date string to be a string composng year and numbe of days
% 
% Usage: output=date2doy('20151212');
%
% by Kang Wang on 01/04/2017
%
yr_str=input(1:4);
month_str=input(5:6);
day_str=input(7:8);

yr=str2num(yr_str);
month=str2num(month_str);
day=str2num(day_str);

lp=isleapyr(yr);

if (lp>0);
    switch month;
        case 1
         iday = day;
        case 2
         iday = 31 + day;
        case 3
          iday= 60 +day;
        case 4
          iday = 91+day;
        case 5
          iday = 121+day;
        case 6
          iday = 152+day;        
        case 7
          iday = 182+day;
        case 8
          iday = 213+day;          
        case 9
          iday= 244+day;
        case 10
          iday = 274+day;
         case 11
          iday = 305+day;
          
         case 12
           iday = 335+day;
           
        otherwise
            error('Wrong input format!')
    end
else

     switch month;
        case 1
         iday = day;
        case 2
         iday = 31 + day;
        case 3
          iday= 59 +day;
        case 4
          iday = 90+day;
          
        case 5
          iday = 120+day;
        case 6
          iday = 151+day;
        case 7
          iday = 181+day;
        case 8
          iday = 212+day;
        case 9
          iday= 243+day;
        case 10
          iday = 273+day;
         case 11
          iday = 304+day;
         case 12
           iday = 334+day;
        otherwise
            error('Wrong input format!')
     end
end
output=[yr_str,sprintf('%03d',iday)];