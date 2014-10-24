%% Version 3.1  5/9/2013
% nullspace seperation of retina
% im_out: original image, high frequencies deleted, masked with a round
% res_out: image in complement of nullspace
% nn_out: image in nullspace
% im: original 2-D image
% slopeRFsz: rate of change of size of RF
% nbands: number of bands in frequency domain
% rfFreqFalloff: width of RF decay region in frequency domain, in octaves
function [im_out,res_out,nn_out] = retina_null(im,slopeRFsz,nbands,rfFreqFalloff)
    
    % debug on-off
    debug = 0;
    dsp = 1; %display bands / results
    
    % check input image
    if size(im,1) ~= size(im,2)
        fprintf('error: width and height of input must be the same.\n');
        return;
    end
    if mod(size(im,1),2) == 1
        fprintf('error: image size must be even number.\n');
        return;
    end
    
    % parameters
    imSz = size(im,1);
    xMax = (imSz-1)/2;

    % spacial domain lattice
    fprintf('set up parameters and functions...\n');
    xvals = (0:(imSz-1))-xMax;
    [x2d,y2d] = meshgrid(xvals);
    % eccentricity from fovea
    r2d = sqrt(x2d.^2+y2d.^2);
    
    % plot original image
    if 0
        figure(1); clf
        plot(xvals,im); rg = get(gca,'Ylim');
        title('original im'); 
        % pause
    end
    
    %slope of RF size w/ eccentricity,
    % slopeRFsz = 1.7/50;  %** parasol, Perry 84
    % slopeRFsz = 1/50;
    % slopeRFsz = 0.4/50; %slope of RF size w/ eccentricity, midget, Perry 84
    %slopeRFsz = 0.22;  % V1
    
    % set up warp function and inverse (image domain and log domain)
    % image domain in pixels, log domain in units of RF spacing
    fovealR = 1/slopeRFsz; %radius, in im pixels, of transition linear->log
    RFsz = (abs(xvals) < 1/slopeRFsz) + (abs(xvals) >= 1/slopeRFsz).* abs(xvals) * slopeRFsz;

    wf = @(r) (r<=fovealR).*r + (r>fovealR).*fovealR.*(log(r./fovealR)+1);
    dwf = @(r) (r<=fovealR).*1.0 + (r>fovealR).*fovealR./r; 

    if (debug)
      x = [-xMax:xMax]; clf;
      subplot(2,1,1);  plot(x,wf(abs(x))); title('warp fn');
      subplot(2,1,2);  plot(x,dwf(abs(x))); title('1/RFsize');
     % pause
    end

    % receptive field filter, in freq domain, relative to foveal sampling
    rfFreqCut = pi;  %**
%     rfFreqFalloff = 1; %*** width of RF decay region, in octaves

    % frequency response of the filter.
    % since the filter is delta function, its DFT is 1 in [-pi,pi]
    % beyond [-pi,pi], it is arbitrary?
    rfFilt = @(f) (abs(f)<rfFreqCut).*1.0 + ...
             ((abs(f)>=rfFreqCut)&(abs(f)<rfFreqCut*2^rfFreqFalloff)).* ...
             cos((pi/2)*log2(abs(f)/rfFreqCut)/rfFreqFalloff).^2;
    
    % plot FT of RF
    if (debug)
      freqs = rfFreqCut * 2^rfFreqFalloff * [-200:200]/200;
       clf; plot(freqs, rfFilt(freqs)); title('RF filter');
      hold on
      plot(-pi*[1 1], [0 1], 'r');
      plot(pi*[1 1], [0 1], 'r');
      hold off
      %pause
    end

%     nbands = 20; %***

    % freq of lowest band, at cutoff, scaled for
    % the largest (outermost) RF sizebandwidth = (pi/minFreq).^(1/(nbands-1));
    % the lowest frequency that is to be thrown away(by the biggest RF)
    % the first frequency that when it's warped, the periphery reach frequency pi
    minFreq = rfFreqCut*dwf(xMax);  
    % compute the bandwidth
    % divide frequency domain into nbands bands, including low-pass part
    bandwidth = (pi/minFreq).^(1/(nbands-1)); 
    % band-pass filter
    % from one octave lower than center frequency to one octave higher than center frequency
    % has value cos(-pi/2) to cos(pi/2), which is 0 to 1 to 0
    bpFilt = @(f) ((abs(f)>1/bandwidth) & (abs(f)<bandwidth)) .* cos((pi/2)*log2(abs(f))./log2(bandwidth));
    % low-pass filter
    % from center frequency to one octave higher than center frequency
    % has value cos(0) to cos(pi/2), which is 1 to 0
    lpFilt = @(f) (abs(f)<=1)*1.0 + ((abs(f)>1) & (abs(f)<bandwidth)) .* cos((pi/2)*log2(abs(f))./log2(bandwidth));
    
    % plot band-pass filter
    if (debug)
        xf = pi * [-xMax:xMax]/xMax;
        clf; plot(xf,bpFilt(xf/(pi/2)));
        %pause
    end

    % % minNullLev = 0.01;  %*** min multiplier for nullspace
    dwiIm = dwf(abs(r2d));

    % frequency domain lattice
    fvals = 2*pi*(0:(imSz-1))/imSz - pi; % [-pi,0),[0 pi)
    [fx2d,fy2d] = meshgrid(fvals);
    freqs = sqrt(fx2d.^2 + fy2d.^2);

    % pre-processing, deleting frequencies higher than pi
    im_f = fft2(im);
    im_f = fftshift(im_f);
    im_f(freqs > pi) = 0;
    im = ifft2(ifftshift(im_f));  
    
    % wavelet analysis
    fprintf('wavelet analysis...\n');
    res=zeros(size(im));  % will contain complement of nullspace
    fmask = zeros(imSz,imSz,nbands+1);
    smask = zeros(imSz,imSz,nbands+1);
    mask = r2d < imSz/2;
    if(dsp), figure(1); clf; end
    fprintf('band %d: ',nbands);
    for n=0:nbands
        fprintf('%d ',n);
        % first ctrFreq is one octave lower than minFreq
        % so low-pass filter goes to 0 at minFreq
        ctrFreq = minFreq*bandwidth^(n-1);
        % filter to extract frequency band
        % fmask will have format of fft of matlab
        % low frequency at beginning, then reach highest frequency
        % then goes back to low frequency
        if (n==0)
            fmask(:,:,n+1) = lpFilt(freqs/ctrFreq);% + lpFilt((sqrt(2)*2*pi-freqs)/ctrFreq);
            fmask(imSz/2+1,imSz/2+1,n+1) = 1;
        else
            fmask(:,:,n+1) = bpFilt(freqs/ctrFreq);% + (n<nbands)*bpFilt((sqrt(2)*2*pi-freqs)/ctrFreq);
            fmask(imSz/2+1,imSz/2+1,n+1)=0; % fmask = fmask ./ sum(fmask(:));

        end
        band = real(ifft2(ifftshift(fmask(:,:,n+1)).*fft2(im)));
        % spatial mask
        % ctrFreq./dwiIm warps the ctrFreq
        % frequency at fovea would be same as ctrFreq
        % frequency at periphery would be higher than ctrFreq
        % at some point it would be higher than pi, the cutoff frequency
        % rfFilt computes what should be the amplitude change at each frequency
        % so far, frequency less than pi has response 1
        % frequency larger than pi to 2*pi has response cos(0)^2 to cos(pi/2)^2
        smask(:,:,n+1) = rfFilt(ctrFreq./dwiIm); % equivalent to ctrFreq.*RFsz
        %mx = max(max(smask(:,:,n+1)),1);
        %smask = 1 ./ (((1-minNullLev)/(mx*minNullLev))*abs(smask) + 1);
        band = band.*smask(:,:,n+1);
        % filter again
        % inverse wavelet transform
        % has anything to do with cumulative fmask not being 1?
        band = real(ifft2(ifftshift(fmask(:,:,n+1)).*fft2(band)));
        res = res + band;
        if (dsp)
            subplot(2,2,1);
            imagesc(xvals,xvals,smask(:,:,n+1)); axis tight square; title('spatial masks');% set(gca,'Xlim',[-xMax,xMax]); set(gca,'Ylim',[0,1]); 
            subplot(2,2,2);
            imagesc(fvals-pi,fvals-pi,fmask(:,:,n+1)); axis tight square; title('freq bands');% set(gca,'Xlim',[-pi,pi]); 
            subplot(2,2,3);
            imagesc(xvals,xvals,band); colormap gray; axis tight square; title('last masked freq band');  % set(gca,'Xlim',[-xMax,xMax]); 
            subplot(2,2,4);
            imagesc(xvals,xvals,res); colormap gray; axis tight square; title('accumulated sum');% set(gca,'Xlim',[-xMax,xMax]); 
    %         pause(0.1);
        end
    end
    fprintf('\n');
    
    im_out = im.*mask;
    res_out = res.*mask;
    nn_out = (im-res).*mask;

    if (dsp)
        h=figure(2);
        set(h,'position',[0 300 1900 500]);
        subplot(1,3,1); imagesc(xvals,xvals,im_out); axis tight square; title('Original im'); % set(gca,'Xlim',xMax*[-1,1],'Ylim',rg); 
        subplot(1,3,2); imagesc(xvals,xvals,res_out); axis tight square; title('Surviving content');
        subplot(1,3,3); imagesc(xvals,xvals,nn_out); axis tight square; title('Nullspace content');

        colormap gray;
    end
    
    
end