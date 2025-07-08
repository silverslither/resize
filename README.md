# resize

ok image resizing implementation.

general-purpose resampling filters ordered from softest (least halo) to sharpest (most halo): triangle, mitnet, catrom, bspline2i, hamming3, bspline3i, hamming4, omoms3i, omoms7i, hamming8, omoms11i.

other filters: nearest (for pixel art), area (for mipmapping), hermite (for somewhat sharp but haloless downsampling), raw bsplines (for blurring).

if you want lanczos, you don't want lanczos (it has terrible effects in both the time and frequency domain). use hamming instead.

if you provide a smooth filter, the example programs will use an edge-directed halo minimization technique to acheive the same amount of halo as your smooth filter, while being just as sharp as the regular filter you pass in.

try it here: https://silverslither.github.io/resize/js/example

TODO: write-up
