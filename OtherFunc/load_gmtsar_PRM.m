function [varargout]=load_PRM(filename,varargin);

[a,b,c]=textread(filename,'%s %s %s');

for i=1:length(varargin);
    string=varargin{i};
%    string=upper(string);
    found=0;
    for j=1:length(a);
        text_tmp = char(a(j));
        if (text_tmp(1)~='#'&strcmp(a(j),string)>0);
            if(length(str2num(char(c(j))))>0);
                varargout{i}= eval(char(c(j)));
            else
                varargout{i}=char(c(j));
            end
                
            found=1;
        end
    end
    if(~found);
        varargout{i}='Variable not found';
    end
end
