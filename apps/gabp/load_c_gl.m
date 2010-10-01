%script for reading the output of the GaBP GraphLab program into matlab
% returns x = inv(A)*b as computed by GaBP
% returns diag = diag(inv(A)) - an approximation to the main diagonal of
% the inverse matrix of A.
% Written by Danny Bickson, CMU

function [ x,diag ] = load_c_gl( fn, n )
  
       
    F = fopen(fn, 'r');
    
    %write matrix size
    x = fread(F,n,'double');
    diag = fread(F, n, 'double'); 
    
    
end

A=rand(3,4);