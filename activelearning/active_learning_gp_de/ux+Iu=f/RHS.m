function f=RHS(x)

    f = 8.*pi.*x.*cos(8.*pi.*x)+sin(8.*pi.*x)+(1/64).*pi.^(-2).*((-8).*pi.*x.*cos(8.*pi.*x)+sin(8.*pi.*x)); 
    
end