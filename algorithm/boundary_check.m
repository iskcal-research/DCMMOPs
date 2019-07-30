function x=boundary_check(x,lower_bound,upper_bound)
% ±ß½ç·´µ¯
length_bound = upper_bound - lower_bound;
x=(x<lower_bound).*(lower_bound+rem((lower_bound-x), length_bound))+(x>=lower_bound).*x;
x=(x>upper_bound).*(upper_bound-rem((x-upper_bound), length_bound))+(x<=upper_bound).*x;

end

