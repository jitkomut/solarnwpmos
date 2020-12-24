
function[] = printtable_head(M,row_head,col_head)
% PRINTTABLE prints the array M in LaTeX table format
% M is an array of size m x n
% The printed table has size of (m+1) x (n+1) because it contains the row
% and column heads
%
% row_head is a cell array of strings of size m x 1
% col_head is a cell array of strings of size (n+1) x 1
% You can modify the format of printing (float, string, integer, etc in the for loop
% USAGE: 
%
% printtableRowColumnhead(M,row_head,col_head)
    
[m,n] = size(M);
printed_col_head = [];
for k=1:n+1
    if k < n+1
        printed_col_head = [printed_col_head '\\bf ' col_head{k} ' & '];
    else
        printed_col_head = [printed_col_head '\\bf ' col_head{k}];
    end
end
printed_col_head = [printed_col_head ' \\\\ \\hline \n'];



fprintf('\n');
fprintf(['\\begin{tabular}{' repmat('c',1,n+1) '} \\hline \n']);
fprintf(printed_col_head);

    for ii=1:m
        value_row = [row_head{ii} ' & '];
            for jj=2:n+1
                if jj < n+1
                    value_row = [value_row [num2str(M(ii,jj-1),'%2.2f'), ' & ']];
                else
                    value_row = [value_row [num2str(M(ii,jj-1),'%2.2f')]];
                end            
            end
            value_row = [value_row ' \\\\ \n'];
            fprintf(value_row);
        
    end
            
fprintf(' \\end{tabular} \n')
end
