%Project 2 Neumann BCs

close all 
clear all
clc
a = 0;
b = 1;
c = 0;
d = 1;
tol = 1e-8;
N = 10;
for p = 1:6
N = 8*p;
h(p) = (b-a)/N;
h1 = h(p)*h(p);
x = a:h(p):b;
y = c:h(p):d;
u2 = zeros(N+1,N+1);
u1 = ones(N+1,N+1);
usol=zeros(N+1,N+1);

for i = 1:N+1
    u2(1,i) = sin(pi*y(i));
    u2(i,1) = 0;
    u2(i,N+1) = 0;
end

k = 0;
while max(max(abs(u2 - u1))) > tol
    k = k+1;
    u1 = u2;
     for i = 2:N+1
        for j = 2:N
             if i == N+1
                u2(i,j) = (2*u1(i-1,j) + u1(i,j-1) + u1(i,j+1)+2*h(p)*exp(1)*sin(pi*y(j))-h1*f(x(i),y(j)))/4;
             elseif i == 2
                 u2(i,j) = (u1(i-1,j) + u1(i+1,j) + u1(i,j-1) + u1(i,j+1)-h1*f(x(i),y(j)))/4;
            elseif i == N
                u2(i,j) = (u1(i-1,j) + u1(i+1,j) + u1(i,j-1) + u1(i,j+1)-h1*f(x(i),y(j)))/4;
            elseif j == 2
                u2(i,j) = (u1(i-1,j) + u1(i+1,j) + u1(i,j-1) + u1(i,j+1)-h1*f(x(i),y(j)))/4;
            elseif j == N
                u2(i,j) = (u1(i-1,j) + u1(i+1,j) + u1(i,j-1) + u1(i,j+1)-h1*f(x(i),y(j)))/4;
            else
           ut = u2;
             u2(i,j) = 4*(ut(i,j-1) + ut(i,j+1) + ut(i+1,j) + ut(i-1,j))/15 -(ut(i-2,j) + ut(i+2,j) + ut(i,j+2) + ut(i,j-2))/60-((3/15)*h1*f(x(i),y(j)));
            
             end
        end
     end
end



for i=1:N+1
    for j=1:N+1
        usol(i,j)=exp(x(i))*sin(pi*y(j));
    end
end


 error(p) = max(max(abs(u2-usol)));
%N = 2*N;
end

for i = 1: p-1
    order(i) = (log(error(i)/error(i+1)))/log(h(i)/h(i+1));
end

 figure(2);
 mesh(x,y,usol);
 title('The exact solution')
 
 figure(1);
 mesh(x,y,u2);
 title('The approximate solution')
order