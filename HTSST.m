function [tfr,t,f] = HTSST(sig,fs,s,method)
%% input:
%sig:Signal
%fs:Sampling frequency
%s:Window function width
%method:Method of time-frequency transform
%       usage: method = struct('type',type,'order',order,'iteration',iter);
%           type:Type of time-frequency transform(Optional options: STFT¡¢TSST)
%           order:Order number (if type is TSST)£¬Optional options: 0 or 1.
%           iter:Iteration number (if type is TSST)
%% output:
%t:Time
%f:Frequency
%tfr:Time-frequency representation
%% 
    if (nargin > 4)
       error('Input has too many parameters.');
    elseif(nargin == 3)
       method = struct('type','STFT');
    elseif(nargin == 2)
       method = struct('type','STFT');s = 0.015;
    elseif(nargin == 1 || nargin == 0 )
        error('Input missing parameters.');
    end
    if ~isstruct(method)
        error('Input parameter error, please check.');
    elseif ~isfield(method,'type')
        error('The input structure should contain type.');
    else
        type = method.type;
        if (~strcmp(type, 'STFT'))&&(~strcmp(type, 'TSST'))
            error('The input struct field "type" should be STFT or TSST.');
        end
        if (strcmp(type, 'TSST'))
            if (~isfield(method,'order')) order = 1;
            else order = method.order; end
            if (order~=1 && order~=2 && order~=3)
                error('The input struct field "order" should be 1 or 2.');
            end
            if (~isfield(method,'iteration')) iter = 1;
            else iter = method.iteration; end
            if (iter~=fix(iter) || iter<=0)
                error('The input struct field "iteration" should be positive integer.');
            end
        end
    end
	%% 
    [xrow,~] = size(sig);
    if (xrow~=1)
        sig = sig';
    end
	%% 
    N = length(sig); 
    dt = 1/fs;
    tao = (0:N-1)*dt;   
    t = (0:N-1)*dt;
    f = ( 0 : round( (N-1)/2 ) )/N * fs;
    L = length(f);
    thre = sqrt(eps);  % threshold value
    %% 
    gt = @(t) s^(-1/2)*pi^(-1/4).*exp(-t.^2/s^2/2);
    Wx = zeros(L, N);
    for ptr = 1:N
        gh = gt(tao-t(ptr));
        xcpsi = fft( conj(gh) .* sig) * dt ;
        Wx(:, ptr) = xcpsi(1:L);
    end
    %% 
    if strcmp(type, 'STFT')
        tfr = Wx;
    else
        %%
        %% 
        gt = @(t) s^(-1/2)*pi^(-1/4).*exp(-t.^2/s^2/2).*(-1j*t);
        dWx = zeros(L, N);
        for ptr = 1:N
            gh = gt(tao-t(ptr));
            xcpsi = fft( conj(gh) .* sig) * dt;
            dWx(:, ptr) = xcpsi(1:L);
        end
        if order == 1
            %%
            GroupDelay =t- real((1i)*dWx./Wx);
        elseif order==2
            %% 
            %% 
            gt = @(t) s^(-1/2)*pi^(-1/4).*exp(-t.^2/s^2/2).*(1j*t/s/s);
            wWx = zeros(L, N);
            for ptr = 1:N
                gh = gt(tao-t(ptr));
                xcpsi = fft(conj(gh) .* sig) * dt;
                wWx(:, ptr) = xcpsi(1:L);
            end
              %% 
         gt = @(t) s^(-1/2)*pi^(-1/4).*exp(-t.^2/s^2/2).*((s.^2-t.^2)/s^4);
            wwWx = zeros(L, N);
            for ptr = 1:N
                gh = gt(tao-t(ptr));
                xcpsi = fft(conj(gh) .* sig) * dt;
                wwWx(:, ptr) = xcpsi(1:L);
            end
            %%            
          Numerator = -Wx.*wWx;
          Denominator = wwWx.*Wx-wWx.*wWx;
            GroupDelay =t- real((1i)*Numerator./Denominator);
        elseif order==3
        %%
            %% 
            gt = @(t) s^(-1/2)*pi^(-1/4).*exp(-t.^2/s^2/2).*(1j*t/s/s);
            %gf = @(w) sqrt(2 * s)*pi^(1/4)*exp(-(s * w).^2/2).*w;
            wWx = zeros(L, N);
            for ptr = 1:N
                gh = gt(tao-t(ptr));
                %gh = gf(tao-t(ptr));
                xcpsi = fft(conj(gh) .* sig) * dt;
                wWx(:, ptr) = xcpsi(1:L);
            end
        %% 
        gt = @(t) s^(-1/2)*pi^(-1/4).*exp(-t.^2/s^2/2).*((s.^2-t.^2)/s^4);
        % gf = @(w) sqrt(2 * s)*pi^(1/4)*exp(-(s * w).^2/2).*w.^2;
            wwWx = zeros(L, N);
            for ptr = 1:N
                gh = gt(tao-t(ptr));
                %gh = gf(tao-t(ptr));
                xcpsi = fft(conj(gh) .* sig) * dt;
                wwWx(:, ptr) = xcpsi(1:L);
            end   
          %% 
        %gt = @(t) s^(-1/2)*pi^(-1/4).*exp(-t.^2/s^2/2).*((1j)*((3*t.*s^2-t.^3)/s^6));
        gt = @(t) s^(-1/2)*pi^(-1/4).*exp(-t.^2/s^2/2).*((1i)*(3*t./s^4-t.^3/s^6));
            wwwWx = zeros(L, N);
            for ptr = 1:N
                gh = gt(tao-t(ptr));
                %gh = gf(tao-t(ptr));
                xcpsi = fft(conj(gh) .* sig) * dt;
                wwwWx(:, ptr) = xcpsi(1:L);
            end
              %% 
        gt = @(t) s^(-1/2)*pi^(-1/4).*exp(-t.^2/s^2/2).*(t.^4/s^8-6*t.^2/s^6+3/s^4);
           wwwwWx = zeros(L, N);
            for ptr = 1:N
                gh = gt(tao-t(ptr));
                %gh = gf(tao-t(ptr));
                xcpsi = fft(conj(gh) .* sig) * dt;
                wwwwWx(:, ptr) = xcpsi(1:L);
            end
         %% 
            Numerator=2*wWx.*wwwWx.*wWx+wwWx.*wwwWx.*Wx-wwWx.*wwWx.*(2*wWx)-wWx.*wwwwWx.*Wx;
            Denominator=Wx.*wwwwWx.*wwWx+wWx.*wwwWx.*wwWx+wWx.*wwwWx.*wwWx-wwWx.*wwWx.*wwWx-wWx.*wwwwWx.*wWx-Wx.*wwwWx.*wwwWx;
            GroupDelay =t- real((1i)*Numerator./Denominator);
          elseif order==4
        end
        
        GroupDelay = 1 + GroupDelay/dt; %    
        GroupDelay = round( round(GroupDelay*2)/2 );
        %% Synchrosqueezing STFT
        tfr = zeros(L,N);
        E = sum(sum(abs(Wx)));
        for prt=1:L
            for b=1:N
                if ( abs(Wx(prt,b)) > thre*E )
                    m = min(max(round(GroupDelay(prt,b)),1),N);
                    tfr(prt, m) = tfr(prt, m) + Wx(prt, b);
                end
            end
        end 
    end
    %% 
	gf = @(w) sqrt(2 * s)*pi^(1/4)*exp(-(s * w).^2/2);
	gh = conj(gf(0));
	tfr = tfr/gh;
end
    
    