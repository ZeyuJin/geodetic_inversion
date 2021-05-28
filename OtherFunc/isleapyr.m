function lp=isleapyr(yr);
% determine if a given year is a leap year
% by Kang Wang on 01/04/2017

if (mod(yr,4)~=0);
    lp=0;
elseif (mod(yr,100)~=0);
        lp=1;
elseif (mod(yr,400)~=0);
            lp=0;
else
            lp=1;
end