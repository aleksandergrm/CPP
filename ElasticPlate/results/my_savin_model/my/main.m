close all;

N = 100;
M = 10*N;
r = linspace(1, 3, N);
fi = linspace(0,2*pi,M);

X = zeros(M,N);
Y = zeros(M,N);
srZ = zeros(M,N);
sfZ = zeros(M,N);
srfZ = zeros(M,N);
for k=1:N
    for l=1:M
        q = r(k)*exp(1i*fi(l));
        [w dw ddw] = omega_zero(q); 
        X(l,k) = real(w);
        Y(l,k) = imag(w);
        [sr sf srf] = stress(q);
        srZ(l,k) = sr;
        sfZ(l,k) = sf;
        srfZ(l,k) = srf;
    end
end

figure
contourf(X,Y,srZ,20)
colorbar
title('Radial stress')
axis equal

figure
contourf(X,Y,sfZ,12)
colorbar
title('Circumferential stress')
axis equal

figure
contourf(X,Y,srfZ,12)
colorbar
title('Shear stress')
axis equal



srV = [];
sfV = [];
srfV = [];

for k=1:N
    q = r(k)*exp(1i*pi/4);
    [sr sf srf] = stress(q);
    srV(end+1) = sr;
    sfV(end+1) = sf;
    srfV(end+1) = srf;
end

figure
plot(r,srfV)


%figure
%plot(r,srV)
%figure
%plot(r,sfV)
%figure
%plot(r,srfV)