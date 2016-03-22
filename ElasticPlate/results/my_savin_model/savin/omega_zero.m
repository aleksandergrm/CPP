function [mp dmp ddmp ] = omega_zero( z )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    mp   = 1/z - 1/6*z^3;
    dmp  = -z^(-2) - 1/2*z^2;
    ddmp = 2*z^(-3) - z;
end

