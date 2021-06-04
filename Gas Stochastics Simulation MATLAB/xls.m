function xls(file,data,tab,rng)

% taskchk=['tasklist /NH /FI "WINDOWTITLE eq Microsoft Excel - ' file '"'];
% [dummy tasks]=dos(taskchk);
% 
% while ~strncmp(tasks,'INFO',4)  
%     pause(2);
%     [dummy tasks]=dos(taskchk);
% end

%warning('off','MATLAB:xls:AddSheet');
try 
    s = 0;
    while s == 0

        warning off MATLAB:xlswrite:AddSheet;
        s=xlswrite(file,data,tab,rng);

    end

catch
        ['Error Writing to ' file ' ' tab]

end
