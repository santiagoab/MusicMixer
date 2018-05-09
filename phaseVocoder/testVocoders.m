%test with pvoc (labrosa)
[waveIn,sr]=wavread('D:\projects\FASTQMUL\silentDisco\sounds\03-FreakFandango-GypsySong.wav', [1000000 2000000]); 
%wavwrite(d,sr,'D:\projects\FASTQMUL\silentDisco\sounds\shortexcerpt.wav');

y=pvoc(waveIn,.75,1024);

wavwrite(waveIn,sr,'D:\matlab\FAST\sounds\Original.wav');
wavwrite(y,sr,'D:\matlab\FAST\sounds\timeStretchedFFTBlock.wav');

%% test synthSignal = pl_phaseVocoder_variable_analysis_hop(signal, tsm_factor, winSamps)
tsm_factor=0.75;
winSamps=2048; %default
synthSignal = pl_phaseVocoder_variable_analysis_hop(waveIn, tsm_factor, winSamps);
wavwrite(synthSignal,sr,'D:\matlab\FAST\sounds\timeStretchedPhaseLockDorran.wav');
