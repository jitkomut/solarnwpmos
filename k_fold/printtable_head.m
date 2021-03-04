
function[] = printtable_head(M,row_head,col_head,printFORMAT,varargin)
% PRINTTABLE prints the array M in LaTeX table format
% M is an array of size m x n
% The printed table has size of (m+1) x (n+1) because it contains the row
% and column heads
%
% row_head is a cell array of strings of size m x 1
% col_head is a cell array of strings of size (n+1) x 1
% 
% printFORMAT is the string of fprintf format, e.g, printFORMAT = '%2.2f' or '%.2e' 
% if printFORMAT of each column is not the same, specify printFORMAT has string cell, e.g, 
%     printFORMAT = {'%2.2f','%.2e','%2.2f','%d'}
% 
% VARARGIN is the binary matrix having the same size as M, if entry is '1' then put \bf (bold face) code in latex
% This matrix should be generated beforehand under the user criterion of interest.
% e.g, when M_{ij} is maximum compared to other entries in the same column
% (or row) we put '1'
% 
% USAGE: 
%
% printtable_head(M,row_head,col_head,'%.2f')
% printtable_head(M,row_head,col_head,'%.2f',BFcode)
% printtable_head(M,row_head,col_head,{'%.2f','%d'}); % suppose M has two
% columns and the values are shown in different printed settings
    

[m,n] = size(M);
if iscell(printFORMAT)
    FORMAT = printFORMAT ;
else
    for jj= 1:n
        FORMAT{jj} = printFORMAT ;
    end
end

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
                    if isempty(varargin)
                        value_row = [value_row [num2str(M(ii,jj-1),FORMAT{jj-1}), ' & ']];
                    elseif varargin{1}(ii,jj-1)
                        value_row = [value_row ['\\bf ' ,num2str(M(ii,jj-1),FORMAT{jj-1}), ' & ']]; % put \bf code
                    else
                        value_row = [value_row [num2str(M(ii,jj-1),FORMAT{jj-1}), ' & ']];
                    end
                else

                    if isempty(varargin)
                        value_row = [value_row [num2str(M(ii,jj-1),FORMAT{jj-1})]];
                    elseif varargin{1}(ii,jj-1)
                        value_row = [value_row ['\\bf ' ,num2str(M(ii,jj-1),FORMAT{jj-1})]];
                    else
                        value_row = [value_row [num2str(M(ii,jj-1),FORMAT{jj-1})]];
                    end
                end            
            end
            value_row = [value_row ' \\\\ \n'];
            fprintf(value_row);
        
    end
            
fprintf(' \\end{tabular} \n')
end
