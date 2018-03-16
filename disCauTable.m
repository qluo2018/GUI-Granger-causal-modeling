function disCauTable(Table, colheads, xLables, TGC_setting, GC_OP)
for i = 1 : length(xLables)
    if length(xLables{i}) < 20
        temp = xLables{i};
        while length(temp)<20
            temp = [temp, ' '];
        end
        rowheads(i,:) = temp;
    end
end
% Create cell array version of table
atab = num2cell(Table);
for i=1:size(atab,1)
   for j=1:size(atab,2)
      if (isinf(atab{i,j}))
         atab{i,j} = [];
      end
   end
end
atab = [cellstr(strjust(rowheads, 'left')), atab];
atab = [cellstr(strjust(colheads, 'left'))'; atab];
digits = [-1 -1 -1 -1 2 4];
if GC_OP == 3
    if TGC_setting.SP_OP==1
        if  TGC_setting.TW_OP==1
            wtitle = 'Spatio-temporal GC';
        elseif TGC_setting.TW_OP==2
            wtitle = 'Spatial GC with fixed time window';
        else
            wtitle = 'Spatial GC';
        end
    else
        if  TGC_setting.TW_OP==1
            wtitle = 'Temporal GC';
        elseif TGC_setting.TW_OP==2
            wtitle = 'Fixed time window GC';
        else
            wtitle = 'classic GC';
        end
    end
elseif GC_OP == 1
    wtitle = 'conditional GC';
elseif GC_OP == 2
    wtitle = 'partial GC';
elseif GC_OP == 4
    wtitle = 'GCSDN';
else
    'error in disCauTable'
end

ttitle = 'Causality Table';
tblfig = statdisptable(atab, wtitle, ttitle, '', digits);
set(tblfig,'tag','table');