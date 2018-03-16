function data = IntegrateDifference(source, d)

data = source;

for i=1:d
    l=size(data,2);
    data = data(:,1:l-1)-data(:,2:l);
end;