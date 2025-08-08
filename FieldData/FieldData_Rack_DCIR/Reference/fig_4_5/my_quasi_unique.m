function qunq_vec = my_quasi_unique( vec )
% MY_QUASI_UNIQUE This function transforms a unique vector into a quasi-unique
% vector, summing up terms of the order of 10^-9 (which must be neglegible for the application)

% Correction
corr = 10^-9;

% Find unique values
[unq_vec,index] = unique(vec,'stable');

% Construction of the quasi-unique vector
k = 1;
for i = 1:length(unq_vec)-1
    for j = 1:index(i+1)-index(i)
        qunq_vec(k) = unq_vec(i) + j*corr;
        k = k+1;
    end
end

% Dealing with the last element
if length(vec)-index(end) > 0
    for j = 1:length(vec)-index(end)
        qunq_vec(k) = unq_vec(end) + j*corr;
        k = k+1;
    end
end
qunq_vec(k) = unq_vec(end);

end