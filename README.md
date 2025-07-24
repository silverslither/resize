# voir

Image resizing implementation in both C99 and JS.

The following are general-purpose resampling filters ordered from softest (least halo) to sharpest (most halo):
- triangle ("linear")
- mitnet ("smooth cubic")
- catrom ("sharp cubic")
- bspline2i (similar to Magic Kernel Sharp 2021)
- hamming3 (similar, but superior to Lanczos3)
- bspline3i
- hamming4
- omoms3i
- omoms7i
- hamming8
- omoms11i

Other filters which are implemented:
- nearest (for pixel art)
- area (for "fair" resampling)
- hermite (for somewhat sharp but haloless downsampling)
- bspline2, bspline3 (for blurring)

Note: You may have noticed that the above lists do not include the Lanczos filters. This is because it exhibits terrible effects in both the time and frequency domain. The Hamming filters are superior windowed-sinc alternatives.

The example programs (a CLI tool, and a simple website) show basic usage of the library. Documentation is available as comments, which should be picked up by your LSP (Doxygen for C, JSDoc for JS).
The example programs also implement an edge-directed halo minimization technique to acheive the same amount of halo as a "smooth" filter, while being just as sharp as the primary filter.

Try the JS example: https://silverslither.github.io/voir/js/example
