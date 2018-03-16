function plotmultiplyTS(data, Ylabel, row, col)
color = colormap(lines);
[T,k] = size(data);
inter = max(max(data));
for i = 1 : k
   subplot(row,col,i)
%    plot(1:T, 2 * inter * (i-1) + data(:,i), 'Color', color(mod(i*10,64)+1,:));
   plot(1:T, data(:,i), 'Color', color(mod(i*10,64)+1,:));
   xlim([0,T]);
%    ylim([min(data(:,i))-0.1, max(data(:,i))+0.1])
%    Ytick(i) = 2 * inter * (i-1);
   if nargin < 2
       Ylabel{i} = ['TS',num2str(i)];
   end
   ylabel(Ylabel{i})
%    set(gca, 'YTick', [] , 'XTick', [])
end
% hold off
xlabel('time')

%%
figure('name', 'ACF')
for i = 1 : k
   subplot(row,col,i)
   autocorr(data(:,i), 40);
   if nargin < 2
       Ylabel{i} = ['TS',num2str(i)];
   end
   ylabel(Ylabel{i})
   set(gca, 'YTick', [] , 'XTick', [])
end
%%
figure('name','ACF of squared sample')
for i = 1 : k
   subplot(row,col,i)
   autocorr(data(:,i).^2, 40);
   if nargin < 2
       Ylabel{i} = ['TS',num2str(i)];
   end
   ylabel(Ylabel{i})
   set(gca, 'YTick', [] , 'XTick', [])
   title('Squared Sample Autocorrelation Function')
end
