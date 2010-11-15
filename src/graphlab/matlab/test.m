function a = test(b) %#eml
%a = sqrt(b)
a.a = double([]);
a.b = int8([1,2,3]);
eml.varsize('a.a', 'a.b');
a.a = chol(b.a);

temp.a = length(b.b(1).c);
temp.b = 20;
a.e = [temp,temp];
a.d = b.b(1).c;
eml.varsize('a.e');
% a.b = int8([]);
% a.c = int32([]);
% a.d = double([]);
% 
% temp.a = 10;
% temp.b = 20;
% a.e = [temp];
% a.f = logical([]);
% a.g = uint16([]);
% a.h = complex(uint32([]));
% a.k = single([]);
% 
% eml.varsize('a.a');
% eml.varsize('a.b');
% eml.varsize('a.c');
% eml.varsize('a.d');
% eml.varsize('a.e');
% eml.varsize('a.f');
% eml.varsize('a.g');
% eml.varsize('a.h');
% eml.varsize('a.k');
% 
% 
% a.d = rand([5,5]);
% 
% 
