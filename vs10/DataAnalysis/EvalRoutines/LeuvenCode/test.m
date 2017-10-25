
tic
x= zeros(1, 1000);

y = zeros(1, 1000);
z = zeros(1, 1000);

for k = 2:1000
   x(k) = x(k-1) + 5;
   y(k) = x(k-1) + 5;
   z(k) = x(k-1) + 5;
end
toc
