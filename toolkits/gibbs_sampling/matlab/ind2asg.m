function asg = ind2asg(siz, ndx)
n = length(siz);
asg = zeros(1,n);
ndx = ndx - 1;
for i = 1:n
  asg(i) = mod(ndx, siz(i));
  ndx = floor(ndx / siz(i));
end
asg = asg + 1;
end
