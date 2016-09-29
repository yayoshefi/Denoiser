function eps = setEpsilon (Type)
global Parameter

if nargin==0;	Type=Parameter.CoOc.Type;   end

switch Type
    case 'CC'
        eps=0.005;
    case 'JP'           % matrix is normalized row & col
        eps=1e-9;
    case 'MI'
        eps=5e-7;     
    case 'PMI'
         eps=0;          %TBD
end

if nargin==0;	Parameter.CoOc.epsilon = eps;   end

end