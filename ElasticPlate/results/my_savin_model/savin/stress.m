function [sr sf srf] = stress( z )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    
    [w dw ddw] = omega_zero(z);
    [phi dphi ddphi] = phi_zero(z);
    [psi dpsi] = psi_zero(z);
    
    mPhi = dphi/dw;
    dmPhi = (ddphi*dw - dphi*ddw)/(dw^2);
    dmPsi = dpsi;
    
    A = 4*real(mPhi);
    B = 2*z^2/(abs(z)^2*conj(dw))*(conj(w)*dmPhi + dmPsi);
    
    sr  = (A - real(B))/2;
    sf  = (A + real(B))/2;
    srf = imag(B);
end

