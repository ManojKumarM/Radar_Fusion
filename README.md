Steps to find XYZ coordinates for a Milli-meter wave RADAR data:

**Radar Used : Infineon BGT60TR13C**

1. Data Collection(100 frames,128 chirps(0.5msec/chirp),64 samples,1 Tx,3 Rx,Total duration = 10sec)
2. Range FFT -> 1D-FFT
3. Velocity FFT -> 2D-FFT
4. CACFAR -> Cell Averaging CFAR
5. Find Azimuth and Elevation Angles based on the Phase Differences between Receiver.
6. XYZ Coordinate Estimation.
7. Plot Animation
8. Plot of X -> Along Azimuth; Y -> Along Range(0 to 120cm); Z -> Along Elevation.
