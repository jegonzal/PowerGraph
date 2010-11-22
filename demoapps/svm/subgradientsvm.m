function [w,b] = subgradientsvm(x,y)
[nd, dim] = size(x);
w=zeros(dim,1);
b = 0;

for j = 1:10000
    stepsize = 1.0/sqrt(j);
    for i = 1:nd
        [dw,db] = wbgradient(x(i,:),y(i),w,b,1.0);
        w = w - stepsize * dw;
        b = b - stepsize * db;
    end
  %  pause    
end

end