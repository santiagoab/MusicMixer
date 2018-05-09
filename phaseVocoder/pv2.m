[d, fs, bps] = wavread ([dir fileName '.wav']);
printf ('Input duration = %.2f s\n', rows (d)/fs);

stretch = 8;
windowsize = round (0.25 * fs);

step = round ((windowsize/2)/stretch);

% original window
fwin = @(x) (1-x.^2).^1.25;
win = fwin (linspace (-1, 1, windowsize));

%win = hanning (windowsize)';

% build index
ind = (bsxfun (@plus, 1:windowsize, (0:step:(rows(d)-windowsize))'))';
cols_ind = columns(ind);

%Only use left channel
left_seg = d(:,1)(ind);
clear d ind;

%Apply window
left_seg = bsxfun (@times, left_seg, win');

%FFT
fft_left_seg = fft (left_seg);
clear left_seg

%keyboard

% overwrite phases with random phases
fft_rand_phase_left = fft_left_seg.*exp(i*2*pi*rand(size(fft_left_seg)));
clear fft_left_seg;

ifft_left = ifft (fft_rand_phase_left);
clear fft_rand_phase_left;

## window again
ifft_left = bsxfun (@times, real(ifft_left), win');

## restore the windowed segments with half windowsize shift
restore_step = floor(windowsize/2);
ind2 = (bsxfun (@plus, 1:windowsize, (0:restore_step:(restore_step*(cols_ind-1)))'))';
left_stretched = sparse (ind2(:), repmat(1:columns (ind2), rows(ind2), 1)(:), real(ifft_left(:)), ind2(end, end), cols_ind);
clear ind2 ifft_left win;

left_stretched = full (sum (left_stretched, 2));

## normalize
left_stretched = 0.8 * left_stretched./max(left_stretched);
printf ("Output duration = %.2f s\n", rows (left_stretched)/fs);
wavwrite (left_stretched, fs, bps, "streched.wav");
system("aplay streched.wav")