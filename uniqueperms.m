function P = uniqueperms(V)
% UNIQUEPERMS - all unique permutations of elements in a vector 
% 
%    P = uniqueperms(V) returns all unique permutations of the N elements
%    in the vector V. P is an array with M rows (see below) and each row of
%    P is unique. The rows are in lexicographic order. V can be a numeric
%    array, or cell array of strings.
%
%    Example:
%       uniqueperms([4 2 1 2]) % returns 12 unique permutions:
%       % [ 1 2 2 4 ; 1 2 4 2 ; 1 4 2 2 ; ... ; 4 2 1 2 ; 4 2 2 1]
%       % perms([3 2 1 2]) will return all 24 permutations of 4 elements.
%       uniqueperms([1 1 1 1 99]) % a 5-by-5 array, rather than a
%                                 % 120-by-5 array with multiple duplications
%
%    This function does not rely on perms to do the job. Similar
%    functionality can be obtained using unique(perms(V),'rows'), but this
%    will create a possibly large intermediate array with N! rows. When V
%    has J unique elements, each of which is repeated Kj times (j=1:J),
%    then the number of unique permutations is N!/(K1!*..*KJ!)
%
%    See also PERMS, NCHOOSEK 
%             NEXTPERM, PERMPOS, PERMSK, PERMN
%
%    Created by Jos (10584)
% tested in Matlab 2019b
% version 2.0 (aug 2019)
% (c) Jos van der Geest
% Matlab File Exchange Author ID: 10584
% email: samelinoa@gmail.com
% History
% 1.0 (aug 2019) created after a File Exchange submission
%     based on an algorithm found at https://en.wikipedia.org/wiki/Permutation
%     and my function NEXTPERM
% 2.0 (aug 2019) optimized code, and implemented two trivial cases
% a straightforward implementation is calling perms with unique on its
% rows, which may produce an unnecessary huge intermediate array. 
% see the commented check below
N = numel(V) ;
P = reshape(V, 1, N) ; % assume trivial permution
if N < 2, return ; end
% convert V to its unique values; work on the indices X into V which allows
% for non-numerical inputs, like cell arrays of chars. Catch the error 
try
    [W, ~, IX] = unique(V) ;
catch ME
    error('UNIQUEPERMS:InputClass', ME.message)
end
if numel(W) == 1, return ; end
% something to do! we use the indices to work on
IX = sort(reshape(IX, 1, N)) ; % we need a sorted row vector
% number of replications of each index
K = histcounts(IX, 1:IX(end)+1) ; 
% number of unique permutations is M = N! / (K1!*...Kn!)
F = cumprod([1 1 2:N]) ; % trick to get all needed factorials (note that all K < N)                                
M = F(N+1) / prod(F(K+1)) ; % number of unique permutations
% pre-allocation using the first permutation. IX is sorted
P = repmat(IX, M, 1) ;
% loop to get the next lexicographic permutation from the previous one
% (see also NEXTPERM on the File Exchange)
Q = P(1, :) ; % previous permution
for k = 2:M
    % find the last element Q(i) that is larger than Q(i-1)
    i = find(Q(2:end) > Q(1:end-1), 1, 'last') ;
    % find the last element Q(j) that is larger than Q(i)
    j = find(Q(i) < Q, 1, 'last') ;
    Q([i j]) = Q([j i]) ;  % swap the two and ...
    Q(i+1:N) = Q(N:-1:i+1) ; % reverse the sequence after the i
    P(k,:) = Q ; % store the new permution, 
    % P0 is the previous one for the next loop
end
% when k = M, all indices are in reversed lexigraphic order: done!
% % CHECK: using unique and perms
% P2 = perms(1:N) ; 
% P2 = unique(IX(P2), 'rows') ;
% isequal(P, P2) 
% replace the indices in P by the element values
P = W(P) ;
end