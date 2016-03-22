function [phi dphi ddphi] = phi_zero( z )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

    phi   = 1/(4*z) + 3/7*z + 1/24*z^3;
    dphi  = -z^(-2)/4 + 3/7 + 1/8*z^2;
    ddphi = 1/2*z^(-3) + 1/4*z;
end

