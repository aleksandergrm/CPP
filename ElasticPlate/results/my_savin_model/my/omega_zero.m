function [mp dmp ddmp ] = omega_zero( z )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    mp   = -1/(6*z^3) + z;
    dmp  = 1 + 1/(2*z^4);
    ddmp = -2/z^5;
end

