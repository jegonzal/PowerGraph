function ndx = asg2ind(siz, asg)
multiple = [1, cumprod(siz(1:end-1))];
assert(isempty(find(asg > siz, 1)));
ndx = sum(multiple .* (asg - 1)) + 1;
end
