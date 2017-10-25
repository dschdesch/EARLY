function [recside,isflipped,wvOUT] = flipchannels(AnNum,contrachan,wvIN)

if contrachan==1
    recside = 'Right';
else
    recside = 'Left';
end

if ~strcmp(AnNum,'L5004') && strcmp(recside,'Right')
    %Flip the channels
    isflipped = 1;
    wvOUT = wvIN(:,[2 1]);
else
    isflipped = 0;
    wvOUT = wvIN;
end
