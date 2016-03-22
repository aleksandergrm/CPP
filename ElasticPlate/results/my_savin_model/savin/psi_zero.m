function [psi dpsi ddpsi] = psi_zero( z )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    psi = -(1/(2*z)) - (91*z - 78*z^3)/(84*(2 + z^4));
    dpsi = 1/(2*z^2) + (z^3*(91*z - 78*z^3))/(21*(2 + z^4)^2) - (91 - 234*z^2)/(84*(2 + z^4));
    %ddpsi = -(1/z^3) - (8 z^6 (91 z - 78 z^3))/(21 (2 + z^4)^3) + (2 z^3 (91 - 234 z^2))/(21 (2 + z^4)^2) + (z^2 (91 z - 78 z^3))/(7 (2 + z^4)^2) + (39 z)/(7 (2 + z^4));
end

