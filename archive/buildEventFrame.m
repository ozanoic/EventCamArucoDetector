function frame = buildEventFrame(events, sensorSize, method)
%BUILDEVENTFRAME Accumulate events into a 2D image frame
%
% frame = buildEventFrame(events, sensorSize, method)
%
% Input:
%   events     - Nx4 matrix [x, y, polarity, timestamp]
%                x, y are 0-indexed; polarity is 0 or 1
%   sensorSize - [height, width] of the sensor
%   method     - 'count'      : event count per pixel (edge map)
%                'integrated' : polarity integration (intensity approx.)
%                'sae'        : surface of active events (latest timestamp)
%
% Output:
%   frame - height x width double matrix

if nargin < 3, method = 'count'; end

height = sensorSize(1);
width  = sensorSize(2);

% Convert 0-indexed to 1-indexed pixel coordinates
x = round(events(:,1)) + 1;
y = round(events(:,2)) + 1;

% Keep only valid coordinates
valid = (x >= 1) & (x <= width) & (y >= 1) & (y <= height);
x = x(valid);
y = y(valid);

idx = sub2ind([height, width], y, x);

switch lower(method)
    case 'count'
        frame = accumarray(idx, 1, [height*width, 1], @sum, 0);

    case 'integrated'
        pol = 2 * events(valid, 3) - 1;   % map {0,1} -> {-1,+1}
        frame = accumarray(idx, pol, [height*width, 1], @sum, 0);

    case 'sae'
        ts = events(valid, 4);
        frame = accumarray(idx, ts, [height*width, 1], @max, 0);

    otherwise
        error('buildEventFrame:badMethod', 'Unknown method: %s', method);
end

frame = reshape(double(frame), [height, width]);
end
