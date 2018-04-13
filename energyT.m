function energyTime=energyT(waveIn, framesize, hop)

% framesize=1024;
% hop=256;
tam=length(waveIn);
nFrames=ceil(tam/hop);
energyTime=zeros(1,nFrames);
waveIn = [waveIn; zeros(framesize,1)];

index=1;
for iframe=1:nFrames   
   grain = waveIn(index:index+framesize);
   media=mean(grain);
   norm=grain-media;
   
   energyFrame=sqrt((sum(norm.^2)/framesize));
   %energyFrame=10*log(sum(grain.^2)/hop);
   %energyFrame=sum(abs(grain/hop));
      
   energyTime(iframe)=energyFrame;
   index=index+hop;
end

%last frame is smaller than framesize
% grain = waveIn(index:tam);
% media=mean(grain);
% norm=grain-media;
%    
% energyFrame=10*log(sum(grain.^2)/hop);
%       
% energyTime(iframe)=energyFrame;
