### Visualization of the LISA pattern functions

`lisajous` = LISA + "Lissajous curves"

This repo provides interactive plots that illustrate the detector response of the Laser Interferometer Space Antenna (LISA), a future space gravitational-wave (GW) detector. The pattern functions give rise to a complex-valued factor as a waveform filtered through the detector response is transformed to the frequency domain. This factor can be represented as a closed contour in the complex plane (a sort of a Lissajous curve). The contour is closed, because the pattern functions are periodic with a period of 1 year (the orbital period of the LISA around the Sun). As such, they can also be decomposed in Fourier series with complex coefficients. `lisajous` also shows these coefficients on a separate plot. 
