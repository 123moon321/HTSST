function awav=wavatten(wavelet,q,timeinterval,samplingrate,pcoeff)
wavlength=length(wavelet);
extl=pcoeff*wavlength;
fwav=fft(wavelet,extl);
hextl=ceil((extl+1)/2);
rhfwav=fwav(1:hextl);
hf=(0:hextl-1)'/(extl*samplingrate);
filter=exp(-pi*timeinterval*hf/q);
arhfwav=rhfwav.*filter;
if mod(extl,2)
    alhfwav=conj(flipud(arhfwav(2:end)));
else
    alhfwav=conj(flipud(arhfwav(2:end-1)));
end
afwav=[arhfwav;alhfwav];
awav=real(ifft(afwav));
awav=awav(1:wavlength);