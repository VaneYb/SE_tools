function out = do_proc_in_freq(y)
    [m,n] = size(y);

    % convert to column vector
    if m < n
        y = y';
        [m,n] = size(y);
    end
    
    % save signal length and the number of channels
    Lx = m;    
    Nch = n;
    clear m n
       
% -------------------------
    
    % define frame length FRMSIZE, the number of overlapped-frame in a 
    % window NOVERLAP, window length WINSIZE
    frmsize = 128;
    noverlap = 2;
    winsize = frmsize * noverlap;
    if (Lx < winsize)
        error('Signal length is too short to process.')
    end
   
    % compute a normalized window
    win = sqrt(hamming(winsize,'periodic'));
    winsigma = sqrt(sum(win.^2) / winsize * noverlap);
    win = win / winsigma;
    win = repmat(win,1,Nch);

    
   % add 0s if signal length is not an integer multiple length of winsize
    Nfrm = ceil((Lx - winsize) / frmsize) + 1;
    L = (Nfrm - 1) * frmsize + winsize;
    out = zeros(L,1);
    y = [y; zeros(L - Lx,Nch)];
    
    % do stft
    nfbands = winsize/2+1;
    Ymat = zeros(nfbands, Nch, Nfrm);
    yfrm = [zeros(frmsize,Nch); y(1:(winsize-frmsize),:)];
    for nt = 1:Nfrm
       yfrm = [yfrm(frmsize+1:end,:); y((winsize-frmsize)+...
           (1+(nt-1)*frmsize:nt*frmsize),:)];
       Yfft = fft(yfrm.*win);
       Ymat(:,:,nt) = Yfft(1:nfbands,:);
    end
    
     
%% Do algorithm
    Xvec = Ymat(:,1,:);
    
    
    
    
    
    
    
    
    
    
    
    
    
%%
    % do istft
    for nt = 1:Nfrm
       X = [ Xvec(:,nt); conj(flip(Xvec(2:end-1,nt)))];
       xfrm = real(ifft(X)).*win(:,1);
       out(1+(nt-1)*frmsize:winsize+frmsize*(nt-1)) = ...
           out(1+(nt-1)*frmsize:winsize+frmsize*(nt-1))+xfrm; 
    end
    out = out(1:Lx);
end