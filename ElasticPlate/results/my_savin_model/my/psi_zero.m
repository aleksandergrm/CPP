function [psi dpsi] = psi_zero( z )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    %psi  = (36*z - 91*z^3 - 84*z^5)/(84 + 168*z^4);
    %dpsi = -((-91 + 77*z^2 - 1300*z^4 + 2128*z^6 + 260*z^8 + 1036*z^10)/( 336*(z + 2*z^5)^2));
    
    psi  = -((91 + 77*z^2 - 130*z^4 + 518*z^6)/(336*(z + 2*z^5)));
    dpsi = -((120 - 247*z^2 + 688*z^4 + 130*z^6 + 168*z^8)/(84*(1 + 2*z^4)^2));
end

