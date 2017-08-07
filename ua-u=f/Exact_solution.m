function f=Exact_solution(x)

global ModelInfo;

a = ModelInfo.alpha;

f = (1/2).*exp(1).^((sqrt(-1)*(-2)).*pi.*x).*(((-1)+((sqrt(-1)*(-2)).* ...
    pi).^a).^(-1).*((sqrt(-1)*(-1))+2.*pi)+exp(1).^((sqrt(-1)*4).*pi.* ...
    x).*((-1)+((sqrt(-1)*2).*pi).^a).^(-1).*(sqrt(-1)+2.*pi));

f = real(f);

end