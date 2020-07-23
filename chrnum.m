function a = chrnum(varargin)

x = varargin{1};

if nargin==1 || strcmp(varargin{2},'hg19') || strcmp(varargin{2},'hg38')
    a = {};
    for i = 1:length(x)
        if x(i)<= 22
            a = [a;num2str(x(i))];
        elseif x(i)==23
            a = [a;'X'];
        elseif x(i)==24
            a = [a;'Y'];
        end
    end
elseif strcmp(varargin{2},'mm10')
    a = {};
    for i = 1:length(x)
        if x(i)<= 19
            a = [a;num2str(x(i))];
        elseif x(i)==20
            a = [a;'X'];
        elseif x(i)==21
            a = [a;'Y'];
        elseif x(i)==22
            a = [a;'M'];
        end
    end
elseif strcmp(varargin{2},'dm6')
    chr_names = [{'2L'} {'2R'} {'3L'} {'3R'} {'4'} {'X'} {'Y'}];
    a = {};
    for i = 1:length(x)
        a = [a;chr_names(x(i))];
    end
end



