%%
vdata.a = [1,2,3];
vdata.b = 'hello';
temp.a = rand(5);
temp.b = uint16(ones(4));
vdata.c = [temp,temp];
vdata.d = [temp,temp];
edata = temp;
%%
compile_update_function({'test_update_function'}, vdata,edata , [getenv('HOME') '/graphlab/graphlabapi/release/src/graphlab'], 'b3');