# Event Camera ArUco Marker Detection

Pure event-based ArUco marker detection using blob analysis and scanline transition decoding. Designed for event camera data (DVS/ESIM), without relying on reconstructed frames.

## Method

1. **Temporal windowing** -- slide through time in 1ms steps, trying multiple lookback windows (5ms to 150ms) to handle markers at different speeds
2. **Event image** -- accumulate events in each window into a binary activity mask using `accumarray`
3. **Blob detection** -- `imfill` + watershed + convex hull to find quadrilateral candidates (`detectQuadBlob.m`)
4. **Perspective unwarp** -- order corners (TL, TR, BR, BL), warp to canonical 160x160 grid via `fitgeotrans`
5. **Scanline transition decoding** -- events appear at cell boundaries, not inside cells. Scan boundary strips vertically and horizontally, flip cell color on threshold crossings, starting from black border
6. **Dictionary lookup** -- match 36-bit code against ARUCO_MIP_36h12 dictionary (250 markers), checking 4 rotations x 2 inversions x 2 flips = 48 candidates per quad

Achieves ~92% detection rate on synthetic ESIM data (95.9% with extended windows).

## Files

| File | Description |
|------|-------------|
| `main.m` | Entry point — configure input file and parameters |
| `detectAruco.m` | Core detection engine: multi-window sliding with optional parallelization |
| `detectQuadBlob.m` | Blob-based quadrilateral detection (imfill + watershed + convex hull) |
| `findQuadCandidates.m` | Region filtering, minimum-area rectangle fitting |
| `loadEvents.m` | Load events from .aedat, .mat, or text files |
| `convertEsimTxt2Mat.m` | Convert ESIM simulator text output to .mat |
| `convertSarmadiBin.m` | Read Sarmadi binary event format |
| `visualizeEvents.m` | Event stream visualization (ON/OFF/combined) |
| `generateMarkerVideo.m` | Generate synthetic marker video for ESIM simulation |
| `generateMarkerTexture.m` | Generate high-res marker texture for ESIM planar renderer |

## Requirements

- MATLAB R2020b or later
- Image Processing Toolbox
- Parallel Computing Toolbox (optional, for `parfor` acceleration)

## Usage

1. Set `matFile` and `sensorSize` in `main.m` to point to your event data
2. Run `main.m`
3. Results are saved to a `.mat` file with per-millisecond detection results across all time windows

## Data Format

Events are expected as an Nx4 matrix: `[x, y, polarity, timestamp]`
- `x, y`: 0-indexed pixel coordinates
- `polarity`: 0 (OFF) or 1 (ON)
- `timestamp`: microseconds

## Synthetic Data Pipeline

To generate synthetic test data using [ESIM](https://github.com/uzh-rpg/rpg_esim):

1. Generate marker texture: `generateMarkerTexture.m`
2. Configure ESIM planar renderer with the texture
3. Run ESIM simulation
4. Convert output: `convertEsimTxt2Mat.m`

## Archive

The `archive/` folder contains experimental approaches that were explored during development (Hough line detection, Sarmadi et al. port, OpenCV-style adaptive thresholding, motion-based detection). These are kept for reference but are not part of the final algorithm.

## Dictionary

Uses the ARUCO_MIP_36h12 dictionary (250 markers, 36-bit codes, minimum Hamming distance 12). The dictionary is embedded directly in `main.m` as a sorted array for parfor compatibility.
