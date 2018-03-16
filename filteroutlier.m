function newdata = filteroutlier(data, alpha)
Mark = zeros(size(data));
for i = 1 : size(data,1)
    stdandardvar = std(data(i,:));
    for j = 1 : size(data,2)
        if abs(data(i,j)) > alpha * stdandardvar
            Mark(i,j) = 1;
        end
    end
end
average = mean(data,2);
newdata = data;
for j = 1 : size(data,2)
    if sum(Mark(:,j)) > 0
        newdata(:,j) = average;
    end
end