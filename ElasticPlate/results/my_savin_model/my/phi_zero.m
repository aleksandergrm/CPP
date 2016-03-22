function [phi dphi ddphi] = phi_zero( z )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

    phi   = 1/(24*z^3) + 3/(7*z) + z/4;
    dphi  = 1/4 - 1/(8*z^4) - 3/(7*z^2);
    ddphi = 1/(2*z^5) + 6/(7*z^3);
end

