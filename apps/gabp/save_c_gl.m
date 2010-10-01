%script for exporting system of linear equations of the type Ax = y to
%graphlab format. 
% Written by Danny Bickson, CMU

function [  ] = save_c_gl( fn, A, y, x )
  
    [sA sA1]= size(A);
    sy = length(y);
    sx = length(x);
    assert(sx == sy);
    assert(sy == sA);
    %assert(sA(1) == sA(2));

    [i1,i2,val] = find((A-diag(diag(A))));
    
    F = fopen(fn, 'w');
    
    %write matrix size
    fwrite(F,sA,'int');
    if (sA == sA1)
        fwrite(F, 0, 'int');
    else
        fwrite(F, sA1,'int'); 
    end
    
   % write y (the observation), x (the solution, if known), diag(A) the
   % variance.
    fwrite(F, y,'double');
    fwrite(F, x, 'double');
    fwrite(F, full(diag(A)), 'double');
    
    %write number of edges
    fwrite(F,length(i1),'int');
    fwrite(F,0, 'int'); % pad with zeros for 64 bit offset
     
    %write all edges
    for i=1:length(i1)
       fwrite(F,i1(i),'int');
       fwrite(F,i2(i),'int');
       fwrite(F,val(i),'double');
    end
     
    fclose(F);
    
    %verify written file header
    F = fopen(fn,'r');
    n = fread(F,1,'int',0,'l');
    assert(n == sA);
    fclose(F);
    
end

