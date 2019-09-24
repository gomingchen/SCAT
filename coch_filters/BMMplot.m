function BMMplot(stimIn, Fs, Fc, FcLabel) 
%function BMMplot(stimIn, Fs, Fc, FcLabel)
% BMMPLOT  Basilar Membrane motion plotter
%
% BMMPLOT(BMM, Fs, Fc, labels)
% Input vector must be in column vector format

switch nargin
    case 4
    case 3
        FcLabel = 1:size(stimIn,2);
    case 2
        FcLabel = 1:size(stimIn,2);
        Fc = 1:size(stimIn,2);
    case 1
        FcLabel = 1:size(stimIn,2);
        Fc = 1:size(stimIn,2);
        Fs = 1;
    otherwise
        error('BMMplot failed, too many parameters!')
end

% initialize memory
stimOut = ones(size(stimIn,1),1) * (1:size(stimIn,2));

% add normalized BM motion to offset vector
stimOut = stimIn./(sqrt(2) * max(max(abs(stimIn)))) + stimOut;

% set time axis
timeAx = linspace(0, length(stimIn)/Fs, length(stimOut)).*1e3;

% plot results
plot(timeAx,stimOut,'k');
set(gca, 'YTick',        FcLabel, ...
         'YTickLabel',   Fc(FcLabel),...
         'TickDir',      'out',...
         'Xlim',         [timeAx(1) timeAx(end)],...
         'Ylim',         [0 size(stimIn,2)+1],...
         'Box',          'off');
axis on
xlabel('Time[ms]');
ylabel('Center Frequency[Hz]');
