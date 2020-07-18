function [pmesh, tmesh, bedges] = GetSimpleDebugMesh()

% Return a very simple triangulation to debug and test functions

simple_pmesh = [0 0 ;
                1 0 ;
                0 1 ;
                1 1 ;];
            
simple_tmesh = [1 2 3;
                2 4 3];
            
simple_bedges = [1    1    2    3 ;
                 2    3    4    4 ;
                100   0    0   101];
             
pmesh  = simple_pmesh;
tmesh  = simple_tmesh;
bedges = simple_bedges';

end

