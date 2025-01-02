# resize

ok image resizing implementation.

general-purpose resampling filters ordered from softest (least halo) to sharpest (most halo): triangle, mitnet, catrom, bspline2i, hamming3, bspline3i, hamming4, omoms3i.

other filters: nearest (for pixel art), area (for mipmapping), hermite (for somewhat sharp but haloless downsampling), raw bsplines (for blurring).

if you want lanczos, you don't want lanczos (it has terrible effects in both the time and frequency domain). use hamming instead.

try it here: https://silverslither.github.io/resize/js/example

TODO: write-up
