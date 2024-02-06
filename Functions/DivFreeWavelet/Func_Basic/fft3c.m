function F = fft3c(f)
% Centered fast Fourier Transform for 3D
    F = fftshift(fftn(ifftshift(f)));
end