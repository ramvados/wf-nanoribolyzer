# DRB004 Barbell filter

The default filter retains forward Ftag matches within the right-end
search window and keeps the sequence before the matched tag and flanks:



    Ftag[fw, *, @right(0..250), <<]

## Dataset validation

Compared with the previous left-plus-right rules on the tested dataset:
- Previous rules retained 143,157 unique reads.
- The right-only rule retained 143,128 unique reads.
- All common reads had identical barcode, flank and cut annotations.
- The 29 excluded reads included 11 with minimap2 hits to the yeast
  RDN37 reference in the segment after the tag and flanks.
- No reference hits were reported for their before-tag segments under
  the same mapping settings; short segments may escape detection.

The right-only workflow completed RNA01 and RNA02 with Dorado 2.1.2:
30 successful tasks, 2 cached tasks, no failures or aborted tasks.
Both HTML reports were produced.

This is a conservative choice supported by this dataset, not a general
validation for every library preparation. Custom filters can be selected
with the barbell_filters parameter.
