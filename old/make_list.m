function bitslist = make_list(L, n1)
    % Initialize cell array to hold binary strings
    bitslist = cell(1, L);
    
    for j1 = 0:(L-1)
        % Convert number to binary string with n1 bits
        bitslist{j1+1} = dec2bin(j1, n1);
    end
end
