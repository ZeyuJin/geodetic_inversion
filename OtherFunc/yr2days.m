function idays=yr2days(year);
%  convert the decimal years to number of days since 01/01/0001
% 
% Usage:  idays=yr2days(year);
% Last Updated by Kang Wang on 01/06/2017 
%
Ndata=length(year);
yr=floor(year);

yr_unique=unique(yr);
nyr=length(yr_unique); %number of unique year
fyr=year-yr; %decimal part of the year

idays=zeros(Ndata,1);
for i=1:nyr;
   this_year=yr_unique(i);
   indx_find=find(yr==this_year);
   fyr_find=fyr(indx_find);
   status=is_leapyear(this_year);
   if (status>0);
     days=365.25; 
   else;
     days=365.25;
   end
   
   doy=fyr_find*days;
   
   date_last_year_end=[num2str(floor(this_year)),'0101'];
%   date_last_year_end=[num2str(floor(this_year-1)),'1231'];
   idays_last_year=datenum(date_last_year_end,'yyyymmdd');
   idays_find=idays_last_year+round(doy);
   idays(indx_find)=idays_find;
 
end

