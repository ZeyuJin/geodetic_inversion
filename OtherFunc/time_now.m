function time_string=time_now();
date_string=date;
clock_string=clock;
hr_now=sprintf('%02d',clock_string(4));
min_now=sprintf('%02d',clock_string(5));
sec_now=sprintf('%02.2f',clock_string(6));

string1=date_string;
string2=[num2str(hr_now),':',num2str(min_now),':',num2str(sec_now)];
time_string=[string1,'-',string2,'-PDT'];
