function [xCoordinates,yCoordinates] = freeHandDraw()
title('Draw the energy curve')
hFH = imfreehand('Closed',false);

% Get the xy coordinates of drawing.
xy = hFH.getPosition;

% get rid of imfreehand remnant.
%delete(hFH);
% Overlay what they drew onto the image.

%hold on; % Keep image, and direction of y axis.

xCoordinates = xy(:, 1);
yCoordinates = xy(:, 2);

%plot(xCoordinates, yCoordinates, 'LineWidth', 2, 'MarkerSize', 10);

end