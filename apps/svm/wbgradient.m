function [dw,db] = wbgradient(x,y,w,b,C)
    alpha = y*(x * w - b);
    if (alpha < 1) 
        alpha = -1;
    elseif (alpha >= 1) 
        alpha = 0;
    end
    %alpha
    dw = w + C*y*alpha*x';
    db = - C*y*alpha;
end