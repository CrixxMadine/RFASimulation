function res = QuadratureTriangle2D(f, intyp, a, b, c)

% Numeric integration of normalized triangles in this form 
% [(a,a), (b,a), (a,c)] 

%  c
%  | \
%  |   \
%  |     \
%  |       \
%  |_ _ _ _ _\
%  a          b  


% returns
% res    :=  numerical result of quadrature

% Input paramaters
% f     := function handle you want to integrate or discrete values
% intyp := type of quadratic integration
% -> 0  =  using matlab function
% -> 1  =  triangular quadrature on 3 points 
% -> 2  =  triangular quadrature on 7 points
%  Formula see ref. book from Jung, tabular 4.14


if (intyp == 1)     % quadrature on 3 points

    m1 = (a(1)+b(1))/2;   % middle point bottom edge
    m2 = (a(2)+c(2))/2;   % middle point left edge

    res = (1/6) * ( f(a(1),m1) + f(m2,a(2)) + f(m2,m1) );     

elseif (intyp == 2)  % quadrature on 7 points

    m1 = (a(1)+b(1))/2;   % middle point lower edge
    m2 = (a(2)+c(2))/2;   % middle point left edge

    d1 = (a(1)+b(1))/3;   % first third lower edge
    d2 = (a(2)+c(2))/3;   % first third left edge

    res = (1/120) * ( 3*f(a(1),a(2)) + 3*f(b(1),b(2)) + 3*f(c(1),c(2)) + ...
                      8*f(a(1),m1) + 8*f(m2,a(2)) + 8*f(m2,m1) + 27*f(d1,d2));


else  

    error('Your integration type is not defined');

end  

        
end

