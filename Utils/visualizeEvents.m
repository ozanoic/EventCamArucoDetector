%% Visualize events in 10ms packets
%  Accumulates events per packet and shows ON/OFF/combined images.

clear; close all; clc;

%% ---- Settings ----
matFile = '../Data/Synthetic/MovingCam/moving_events/moving_events.mat';
packetDt = 10000;   % 10 ms in microseconds
sensorSize = [240, 320];  % [H, W]

%% ---- Load ----
fprintf('Loading %s...\n', matFile);
tmp = load(matFile, 'events');
events = tmp.events;
numEvents = size(events, 1);
fprintf('Loaded %d events\n', numEvents);

H = sensorSize(1);
W = sensorSize(2);

%% ---- Group into packets ----
tAll = events(:, 4);
tMin = min(tAll);
tMax = max(tAll);
packetEdges = tMin : packetDt : (tMax + packetDt);
numPackets = length(packetEdges) - 1;
[~, ~, packetIdx] = histcounts(tAll, packetEdges);

fprintf('Packets: %d  (%.0f us each)\n', numPackets, packetDt);

%% ---- Visualize ----
hFig = figure('Name', 'Event Visualization', 'Position', [50 50 1200 400]);

for p = 1:numPackets
    idx = (packetIdx == p);
    if sum(idx) < 1, continue; end

    pktEvents = events(idx, :);
    px  = pktEvents(:,1) + 1;   % 1-indexed
    py  = pktEvents(:,2) + 1;
    pol = pktEvents(:,3);        % 0=OFF, 1=ON

    % Build images
    onImg  = zeros(H, W);
    offImg = zeros(H, W);
    for e = 1:size(pktEvents, 1)
        r = py(e); c = px(e);
        if r < 1 || r > H || c < 1 || c > W, continue; end
        if pol(e) == 1
            onImg(r, c) = onImg(r, c) + 1;
        else
            offImg(r, c) = offImg(r, c) + 1;
        end
    end

    % Combined RGB: ON=red, OFF=blue
    combRGB = zeros(H, W, 3);
    combRGB(:,:,1) = min(onImg / max(onImg(:) + eps), 1);
    combRGB(:,:,3) = min(offImg / max(offImg(:) + eps), 1);

    figure(hFig); clf;

    % subplot(1,3,1);
    % imagesc(onImg); colormap(gca, 'hot'); axis image; colorbar;
    % title(sprintf('ON events (%d)', sum(pol == 1)));
    % 
    % subplot(1,3,2);
    % imagesc(offImg); colormap(gca, 'hot'); axis image; colorbar;
    % title(sprintf('OFF events (%d)', sum(pol == 0)));
    % 
    % subplot(1,3,3);
    imshow(combRGB); axis image;
    title('Combined (R=ON, B=OFF)');

    tStart = packetEdges(p);
    sgtitle(sprintf('Packet %d/%d  |  t=%.3f s  |  %d events', ...
        p, numPackets, tStart/1e6, sum(idx)));
    drawnow;
end

fprintf('Done.\n');
