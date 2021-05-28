function output=doy2date(input);
% convert the  string composng year and numbe of days to string
% 
% Usage: output=doy2date('2015035');
% 
% by Kang Wang on 01/04/2017
%
yr_str=input(1:4);
doy_str=input(5:7);
yr=str2num(yr_str);
doy=str2num(doy_str);
lp=isleapyr(yr);
day_of_year=zeros(13,1);

if (lp>0);
    day_of_yr(1)=0;
    day_of_yr(2)=31;
    day_of_yr(3)=60;
    day_of_yr(4)=91;
    day_of_yr(5)=121;
    day_of_yr(6)=152;
    day_of_yr(7)=182;
    day_of_yr(8)=213;
    day_of_yr(9)=244;
    day_of_yr(10)=274;
    day_of_yr(11)=305;
    day_of_yr(12)=335;
    day_of_yr(13)=366;
else
    day_of_yr(1)=0;
    day_of_yr(2)=31;
    day_of_yr(3)=59;
    day_of_yr(4)=90;
    day_of_yr(5)=120;
    day_of_yr(6)=151;
    day_of_yr(7)=181;
    day_of_yr(8)=212;
    day_of_yr(9)=243;
    day_of_yr(10)=273;
    day_of_yr(11)=304;
    day_of_yr(12)=334;
    day_of_yr(13)=365;
end

if (doy>day_of_yr(1)&doy<=day_of_yr(2));
    imonth=1;
    iday=doy;
elseif (doy>day_of_yr(2)&doy<=day_of_yr(3));
    imonth=2;
    iday=doy-day_of_yr(2);
    
elseif (doy>day_of_yr(3)&doy<=day_of_yr(4));
    imonth=3;
    iday=doy-day_of_yr(3);
    
elseif (doy>day_of_yr(4)&doy<=day_of_yr(5));
    imonth=4;
    iday=doy-day_of_yr(4);
    
elseif (doy>day_of_yr(5)&doy<=day_of_yr(6));
    imonth=5;
    iday=doy-day_of_yr(5);  
    
elseif (doy>day_of_yr(6)&doy<=day_of_yr(7));
    imonth=6;
    iday=doy-day_of_yr(6); 
    
elseif (doy>day_of_yr(7)&doy<=day_of_yr(8));
    imonth=7;
    iday=doy-day_of_yr(7);
    
elseif (doy>day_of_yr(8)&doy<=day_of_yr(9));
    imonth=8;
    iday=doy-day_of_yr(8);
    
elseif (doy>day_of_yr(9)&doy<=day_of_yr(10));
    imonth=9;
    iday=doy-day_of_yr(9);    
    
elseif (doy>day_of_yr(10)&doy<=day_of_yr(11));
    imonth=10;
    iday=doy-day_of_yr(10);
    
elseif (doy>day_of_yr(11)&doy<=day_of_yr(12));
    imonth=11;
    iday=doy-day_of_yr(11);  
    
 elseif (doy>day_of_yr(12)&doy<=day_of_yr(13));
    imonth=12;
    iday=doy-day_of_yr(12);
end

day_str=sprintf('%02d',iday);
mon_str=sprintf('%02d',imonth);
output=[yr_str,mon_str,day_str];
    
    