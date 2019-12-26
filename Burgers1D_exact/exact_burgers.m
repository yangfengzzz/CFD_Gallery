function [ u, status] = exact_burgers(x, t)
%EXACT_BURGERS 

xi1 = fzero(@(xx) 0.3+0.7*sin(xx*pi)+(xx-x)/t, -5);
xi2 = fzero(@(xx) 0.3+0.7*sin(xx*pi)+(xx-x)/t, 5);
if abs(xi1 - xi2)
    xi = xi1;
    status = 0;
else
   val1 = 0.3*xi1 - 0.7*cos(xi1*pi)/pi + 0.7/pi + (x-xi1)^2/2/t;
   val2 = 0.3*xi2 - 0.7*cos(xi2*pi)/pi + 0.7/pi + (x-xi2)^2/2/t;
   if val1 < val2
       xi = xi1;
   else
       xi = xi2;
   end
   status = 1;
end

u = (x-xi)/t;
end

