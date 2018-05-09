
spectrogram(features(1).wave);

spectrum = pwelch(features(1).wave);
plot(10*log10(spectrum))
drawnow;
sound(features(1).wave, features(1).sr);