function res = QuadratureTriangle2D(f, ftype, intyp, a, b, c)

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
% ftype := indicator, whether f is function handle or discrete values
% -> 1 = f is function handle {z = f(x,y)}
% -> 2 = if is vector with discrete values [z_a; z_b; z_c)
% [a b] := left vertical edge 
% [a c] := lower horizontal edge
% 
% intyp := type of quadratic integration
% -> 0  =  using matlab function
% -> 1  =  triangular quadrature on 3 points 
% -> 2  =  triangular quadrature on 7 points
%  Formula see ref. book from Jung, tabular 4.14


if (ftype == 1) % f is a function handle 

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


    else  % matlab library function  -> 
          % actually very slow and inaccurate

        res = integral2(f, a(1), b(1), a(2), c(2));
    
    end  % intyp == ?

elseif (ftype == 2) % f has discrete values on the corners
   
    % Area of triangle with given coordinates 
    % See Wikipedia ...    
    A = [  1   ,  1   ,  1   ;
          a(1) , b(2) , c(2) ; 
          a(2) , b(2) , c(2) ; ];     
    area = 0.5 * abs(det(A));
    
    % Ref: https://uk.mathworks.com/matlabcentral/answers/498300
    %      -how-to-integrate-discrete-values-over-a-surface#comment_781700
    res = 0.5 * (f(1) + f(2) + f(3)) * area;
        
end

