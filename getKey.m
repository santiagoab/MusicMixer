function [key, mode] = getKey(dir)

addpath('C:\Users\SANTI\Desktop\MusicMixer\MIRtoolbox1.7\MIRToolbox')
addpath('C:\Users\SANTI\Desktop\MusicMixer\MIRtoolbox1.7\MIRToolboxDemos')

c = mirchromagram(dir);
ks = mirkeystrength(c);
[k kc ks] = mirkey(ks);
[d,d2] = mirgetdata(k);

key = d;
mode = d2;

% % cambiar 1-12 por notas
% 
% switch d
%     case 1
%         key = 'C';
%     case 2
%         key = 'C#';
%     case 3
%         key = 'D';
%     case 4
%         key = 'D#';
%     case 5
%         key = 'E';
%     case 6
%         key = 'F';
%     case 7
%         key = 'F#';
%     case 8
%         key = 'G';
%     case 9
%         key = 'G#';
%     case 10
%         key = 'A';
%     case 11
%         key = 'A#';
%     case 12
%         key = 'B';
% end
% 
% % cambiar 1-2 por Maj, Min
% switch d2
%     case 1
%         mode = 'Maj';
%     case 2
%         mode = 'Min';
% end

end