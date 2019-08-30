function x = opLadj_1D(y, opt)

l=length(y);

if strcmp(opt.filter,'gradient');
    if strcmp(opt.computation,'fourier')
        h = [1/2, -1/2];
        H = psf2otf(h,size(y));
        x = real(ifft(fft(y).*conj(H)));
    elseif strcmp(opt.computation,'direct')
        x = -([y(1:1)/2,y(2:l-1)/2- y(1:l-2)/2,-y(l-2+1:l-1)/2]);
    end
    
elseif strcmp(opt.filter,'laplacian');
    if strcmp(opt.computation,'fourier')
        h = [1/4, -1/2, 1/4];
        H = psf2otf(h,size(y));
        x = real(ifft(fft(y).*conj(H)));   
    elseif strcmp(opt.computation,'direct')
        x = ([y(1:1)/4,-y(1,1)/2 + y(1,2)/4, y(3:l-2)/4 - y(2:l-3)/2 + y(1:l-4)/4, y(l-3)/4 - y(l-2)/2, y(l-2)/4]);
    end
end