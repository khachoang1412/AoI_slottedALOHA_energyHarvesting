function S = integer_partitions(n,count)
%
% This work is licensed under a 
% Creative Commons Attribution 4.0 International License. 
% https://creativecommons.org/licenses/by/4.0/
%
% Coded by IGOR S. PERETTA - iperetta@gmail.com (2015)
%
% Arguments/inputs:
% - n: the desired nonnegative integer number to be partitioned
% - count: the maximum number of integers to build the desired number
%
% Output:
% - S: a matrix with "count" columns and as many lines as necessary to 
%      provide every possible partition of a nonnegative integer "n" 
%      considering a sum of "count" integer numbers. Related permutations 
%      are not considered here (try "help perms" for that).
%
% Examples of usage:
%         >> integer_partitions(5)
%         
%         ans =
%         
%              5     0     0     0     0
%              4     1     0     0     0
%              3     2     0     0     0
%              3     1     1     0     0
%              2     2     1     0     0
%              2     1     1     1     0
%              1     1     1     1     1
%         
%         >> integer_partitions(5,3)
%        
%         ans =
%         
%              5     0     0
%              4     1     0
%              3     2     0
%              3     1     1
%              2     2     1
%        
%         >> integer_partitions(3,6)
%         
%         ans =
%         
%              3     0     0     0     0     0
%              2     1     0     0     0     0
%              1     1     1     0     0     0
%
% Adapted from "Algorithm ZS1" in
% ANTOINE ZOGHBI and IVAN STOJMENOVIC (1998), "FAST ALGORITHMS FOR 
% GENERATING INTEGER PARTITIONS", International Journal of Computer 
% Mathematics, Volume 70, Issue 2, pages 319-332
% DOI: 10.1080/00207169808804755
%
if nargin == 1
    count = n;
end
if n < 0 || n ~= round(n)
    error('Only nonnegative integers allowed!');
elseif n == 0
    if count == 0
        S = 0;
    else
        S = zeros(1,count);
    end
else
    x = ones(1,n);
    x(1) = n;
    m = 1;
    h = 1;
    M = [x(1:m) zeros(1,n-m)];
    while x(1) ~= 1
        if x(h) == 2 
           m = m + 1;
           x(h) = 1;
           h = h - 1;
        else
           r = x(h) - 1;
           t = m - h + 1;
           x(h) = r;
           while t >= r
               h = h + 1;
               x(h) = r;
               t = t - r;
           end
           if t == 0
               m = h;
           else
               m = h + 1;
               if t > 1
                   h = h + 1;
                   x(h) = t;
               end
           end
        end
        M = cat(1,M,[x(1:m) zeros(1,n-m)]);
    end
    if count > n
        M = cat(2,M,zeros(size(M,1),count-n));
    end
    S = [];
    for i = 1:size(M,1)
        if(sum(M(i,1:count)) == n)
            S = cat(1,S,M(i,1:count));
        end
    end
end