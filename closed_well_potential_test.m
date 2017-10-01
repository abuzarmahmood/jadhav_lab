[X Y] = meshgrid(-100:0.1:100,-100:0.1:100);

all_ones = ones(size(X,1),size(X,2));
Z = zeros(size(X,1),size(X,2));

for x = 0:50
    Z1 =  all_ones(X==x);
    Z1 = (1+exp(-(Z.^2)).^(-1));
    Z = Z+Z1;
    print(x)
end
%%
% Make grid spacing low, the actual potential can have high spacing
% but gradient has low spacing
f = @(x,y,x1,y1,w) (1+exp(-(w.*(x-x1).^2+w.*(y-y1).^2))).^-1;
[x y]=meshgrid(-100:0.1:100,-100:0.1:100);
w = 1/20; %width parameter

x_s = [-50,50];
y_s = [-50,50];
test_line_x = linspace(x_s(1),x_s(2),100);
test_line_y = linspace(y_s(1),y_s(2),100);

% Don't do this...takes too long
z_f = zeros(size(x,1),size(x,2));

for i = 1:length(test_line_x)
    z =feval(f,x,y,test_line_x(i),test_line_y(i),w);
    z_f = z_f + z;
    disp(i)
    
end
%{
z_f2 = zeros(size(x,1),size(x,2));
for y1 = -50:1:50
    z =feval(f,x,y,0,y1);
    z_f = z_f + z;
end


z = z_f + z_f2;
%}

[px py] = gradient(z_f,0.1);
quiver(x,y,px,py)

%%
f = @(x,y,x1,y1) (1+exp(-((x-x1).^2+(y-y1).^2))).^-1;
x = -100:0.1:100;
y = -100:0.1:100;
y1=0;
fun = @(x1) f(x,y,x1,y1);
z_f = feval(fun,[-50,50]);