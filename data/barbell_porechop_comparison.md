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
