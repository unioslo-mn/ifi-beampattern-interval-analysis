function [PSLL, PSLL_idx] = PSLL(Upper, Nominal, axis)
    [ ~, idx ] = max( Nominal );
    min_indx = find( islocalmin(Nominal) == true );
    
    if idx <= min_indx(1)
        start = 1;
        stop = min_indx(1);
    elseif idx >= min_indx(end)
        start = min_indx(end);
        stop = length(axis);
    else % somewhere in between
        for i = 1:length(min_indx)
            if idx <= min_indx(i)
                start = min_indx(i-1); % won't try to access (0) because of first if clause
                stop = min_indx(i); 
                break
            end
        end
    end
    
    Upper(start:stop) = 0;
    [PSLL, PSLL_idx] = max(Upper);
                
end

