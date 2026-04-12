function events = convertSarmadiBin(binFilePath, matFilePath)
%CONVERTSARMADIBIN Convert Sarmadi et al. packets.bin to our .mat format.
%
%  The .bin format (from event-aruco C++ code, Linux 64-bit) is:
%    size_t    num_packets            (8 bytes, uint64)
%    For each packet:
%      int64   timestamp              (8 bytes, packet-level timestamp ms)
%      size_t  num_events             (8 bytes, uint64)
%      int16   packet_source          (2 bytes)
%      int32   packet_ts_overflow     (4 bytes)
%      For each event:
%        uint16  x                    (2 bytes)
%        uint16  y                    (2 bytes)
%        int32   ts                   (4 bytes, event timestamp in us)
%        uint8   polarity             (1 byte, bool: 0=OFF, 1=ON)
%        uint8   validity             (1 byte, bool)
%
%  Output .mat format (compatible with loadEvents.m):
%    events = [x, y, polarity, timestamp]
%      x, y       : 0-indexed pixel coordinates
%      polarity   : -1 (OFF) or +1 (ON)
%      timestamp  : microseconds (int32 from the event)
%
%  Usage:
%    events = convertSarmadiBin('path/to/packets.bin');
%    events = convertSarmadiBin('path/to/packets.bin', 'path/to/output.mat');
%
%  The DVS128 sensor used by Sarmadi is 128×128 pixels.

if nargin < 2, matFilePath = ''; end

fprintf('Reading: %s\n', binFilePath);
fid = fopen(binFilePath, 'rb', 'l');   % little-endian
if fid < 0
    error('Cannot open file: %s', binFilePath);
end

%% ---- Read header -----------------------------------------------------------
num_packets = fread(fid, 1, 'uint64');
fprintf('Number of packets: %d\n', num_packets);

%% ---- Read all packets ------------------------------------------------------
% Pre-scan to estimate total events (read packet headers)
% We'll collect all events into a growing cell array, then concatenate.
allEvents = cell(num_packets, 1);
totalEvents = 0;

for p = 1:num_packets
    % Packet header
    pkt_timestamp    = fread(fid, 1, 'int64');     % ms timestamp
    num_events       = fread(fid, 1, 'uint64');
    packet_source    = fread(fid, 1, 'int16');
    ts_overflow      = fread(fid, 1, 'int32');

    if isempty(num_events) || num_events == 0
        continue;
    end

    % Read all events in this packet at once (10 bytes per event)
    % Layout per event: uint16 x, uint16 y, int32 ts, uint8 pol, uint8 val
    % = 2 + 2 + 4 + 1 + 1 = 10 bytes
    raw = fread(fid, [10, num_events], '*uint8');  % 10 x num_events
    if size(raw, 2) < num_events
        warning('Packet %d: expected %d events but read %d. File may be truncated.', ...
            p, num_events, size(raw, 2));
        num_events = size(raw, 2);
        if num_events == 0, continue; end
    end

    % Parse fields from raw bytes (little-endian)
    x   = double(typecast(reshape(raw(1:2, :), [], 1), 'uint16'));
    y   = double(typecast(reshape(raw(3:4, :), [], 1), 'uint16'));
    ts  = double(typecast(reshape(raw(5:8, :), [], 1), 'int32'));
    pol = double(raw(9, :)');      % 0=OFF, 1=ON
    val = double(raw(10, :)');     % validity

    % Keep only valid events
    validIdx = (val > 0);
    x   = x(validIdx);
    y   = y(validIdx);
    ts  = ts(validIdx);
    pol = pol(validIdx);

    % Convert polarity: 0 → -1, 1 → +1
    pol(pol == 0) = -1;

    nValid = sum(validIdx);
    if nValid > 0
        allEvents{p} = [x, y, pol, ts];
        totalEvents = totalEvents + nValid;
    end

    % Progress
    if mod(p, 500) == 0 || p == num_packets
        fprintf('  Packet %d / %d  (events so far: %d)\n', p, num_packets, totalEvents);
    end
end

fclose(fid);

%% ---- Concatenate and sort by timestamp -------------------------------------
allEvents = allEvents(~cellfun('isempty', allEvents));
events = vertcat(allEvents{:});

% Sort by timestamp
events = sortrows(events, 4);

fprintf('\nTotal valid events: %d\n', size(events, 1));
fprintf('  x range: [%d, %d]\n', min(events(:,1)), max(events(:,1)));
fprintf('  y range: [%d, %d]\n', min(events(:,2)), max(events(:,2)));
fprintf('  t range: [%d, %d] us  (duration: %.1f ms)\n', ...
    min(events(:,4)), max(events(:,4)), ...
    (max(events(:,4)) - min(events(:,4))) / 1000);
fprintf('  polarity: %d OFF, %d ON\n', ...
    sum(events(:,3) == -1), sum(events(:,3) == 1));

%% ---- Save to .mat ----------------------------------------------------------
if ~isempty(matFilePath)
    save(matFilePath, 'events', '-v7.3');
    fprintf('Saved to: %s\n', matFilePath);
end

end
