# Porechop and Barbell comparison

## Scope

This comparison uses the older preprocessing outputs from 2 October 2026,
not the later Dorado 2.1.2 right-only filter test.
Both preprocessing runs used the same untrimmed FASTQ input
(1,107,032 reads). Input SHA256:
63812e9156ee00b73a0546065b29680bb1236dd92a0045f1d61e87e29dd46b50

The Barbell output used the previous left-plus-right Ftag rules.
The current right-only default was tested separately.

Reads were aligned to references/RDN37-1.fa with minimap2 2.28-r1209:
minimap2 -x map-ont -c --secondary=no -t 2

## Overall read yield

| Metric | Porechop | Barbell |
|---|---:|---:|
| Retained reads | 1,031,584 | 142,187 |
| Reads with an rRNA alignment | 57,615 | 6,560 |
| Percentage with an rRNA alignment | 5.59% | 4.61% |

There were 124,722 shared read IDs.
Of the exclusive read IDs, 51,058 Porechop-only reads and 7 Barbell-only
reads had an rRNA alignment.

## Shared reads

| Metric | Porechop | Barbell |
|---|---:|---:|
| Reads with an rRNA alignment | 6,557 | 6,553 |
| Aligned query bases, overlapping intervals counted once | 2,473,762 | 2,474,718 |
| Median aligned query bases per mapped read | 148 | 148 |

Among the 6,553 reads mapped by both methods:
- 6,333 had equal aligned query lengths.
- 220 had greater aligned query lengths after Barbell.
- None had greater aligned query lengths after Porechop.

The four Porechop-only alignments among shared reads belonged to the
previously examined Ftag exception reads. Their rRNA alignments were
after the barcode, whereas the old before-tag cut retained the prefix.

## Where reads were excluded

All 51,058 Porechop-only rRNA-aligned reads lacked an entry in the
matching Barbell anno.tsv. They were not excluded by the downstream
filter rules. Absence from the annotation table does not establish
that a biological barcode was absent.

A further 970 IDs were present in filtered.tsv but absent from the
combined Barbell FASTQ. All had read_start_flank = 0, consistent with
a before-tag cut leaving an empty prefix.

## Interpretation and limitations

Barbell largely preserved the detected rRNA alignments among shared
reads, but its configured barcode selection produced substantially
fewer rRNA-aligned reads overall than Porechop.

This comparison does not establish barcode assignment accuracy,
explain why the missing barcodes were not detected, or establish
equivalence on other datasets. Short reads may fail to align with
the chosen mapping settings.

Barbell is an optional demultiplexing route. Porechop remains the
workflow default. The right-only filter change does not resolve
the much larger difference caused at the annotation stage.

## Boundary-penalty experiment

Barbell 0.3.3 was tested on the same untrimmed input with alpha 0.4
(default) and alpha 0.2. Both used the same right-only Ftag filter.
The automatic flank edit cutoff remained 5.

### Diagnostic subset

The subset contained 200 previously recognized rRNA reads and 200
Porechop-aligned rRNA reads without a previous Barbell annotation.

- Default settings assigned barcodes to the 200 recognized controls only.
- Flank cutoff 7 recovered 11 additional barcode candidates.
- Alpha 0.2 recovered 22 additional barcode candidates.
- Existing barcode/orientation assignments were unchanged.
- An independent exact 18-base barcode-core search found candidates in
  104/200 missed reads and 188/200 recognized controls.
- Five inspected missed examples ended within the barcode sequence.

These observations support incomplete barcode/flank sequences as one
cause of missed annotations, but do not establish assignment accuracy.

### Artificial truncation test

Five variants were generated for each of 199 recognized controls;
one control was excluded because its reported boundaries were unsuitable.
Labels were compared with the original Barbell assignment, not an
independent experimental ground truth. Cuts used reported coordinates.

| Variant, 199 reads each | Alpha 0.4: original assignment retained | Alpha 0.2: original assignment retained |
|---|---:|---:|
| Unmodified | 199 | 199 |
| Right flank removed | 0 | 195 |
| Another 10 barcode bases removed | 0 | 7 |
| Another 20 barcode bases removed | 0 | 0 |

No differing barcode/orientation assignments were observed.
The 199 prefixes with the selected tag removed produced no Ftag hits
under either setting. This is a limited control, not a comprehensive
false-positive assessment.

### Full-input trimming and alignment

Both variants were trimmed and aligned with the same minimap2 settings
used above.

| Metric | Alpha 0.4, right-only | Alpha 0.2, right-only |
|---|---:|---:|
| Unique reads emitted | 142,159 | 346,373 |
| Reads aligned to RDN37 | 6,560 | 12,453 |
| Percentage aligned | 4.61% | 3.60% |
| Aligned query bases | 2,477,354 | 4,317,318 |
| Median aligned query bases per mapped read | 148 | 148 |

There were 5,937 shared aligned reads with identical aligned query
lengths, 623 standard-only aligned reads and 6,516 alpha-0.2-only
aligned reads.

At the filter stage, alpha 0.2 excluded 6,624 previously selected IDs
despite all retaining an Ftag annotation. Of these, 6,545 also had
Fflank annotations; composite patterns are not accepted by the current
single-Ftag rule. The other cases require separate inspection.

Alpha 0.2 increases absolute rRNA yield but also changes pattern
selection. Correct assignment of the additional reads has not been
independently validated. The workflow retains the Barbell default
alpha 0.4; alpha 0.2 remains an experimental setting.
