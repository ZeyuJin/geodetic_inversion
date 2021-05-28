function yr_out = days2yr(days);
% convert the number of days since 01/01/0000 to dicimal yr
% yr=days2yr(days);
% Last update by Kang Wang on 01/06/2017
Ndays=length(days);
yr_out=zeros(Ndays,1);
for i=1:Ndays;
  str_tmp=datestr(days(i),'yyyymmdd');
  yr_str=str_tmp(1:4);
  yr=str2num(yr_str);
  days_till_last_yr=datenum([num2str(floor(yr)),'0101'],'yyyymmdd');
  doy=days(i)-days_till_last_yr;
  
  lp=is_leapyear(yr);
  if (lp>0);
   days_year=365.25;
  else;
   days_year=365.25;
  end

  fyr=doy/days_year;
  yr_out(i)=yr+fyr;
end


