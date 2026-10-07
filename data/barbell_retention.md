# Barbell read retention

Barbell uses its default alpha of 0.4 and the right-end Ftag filter.
Reads present in Barbell's trimmed output retain their trimmed sequences.
All remaining input reads are retained unchanged in a separate
unassigned group. This group does not imply a barcode assignment.

With demultiplexing enabled, unassigned reads are analyzed separately.
The unassigned group is retained even when specific barcodes are selected.
Other barcode groups remain subject to that selection.
Without demultiplexing, trimmed and fallback reads are combined.

## Dataset validation

Run: nauseous_pauling, completed on 2026-10-06 on MOGON.
Dorado 2.1.2 and Barbell 0.3.3 were used.

- Raw input: 1,107,032 reads.
- Barbell trimmed output: 142,159 reads.
- Raw fallback: 964,873 reads.
- Accepted by the filter but absent from trimmed output: 969 reads.
- Not retained by the filter: 963,904 reads.
- Workflow: 46 successful tasks, 2 cached tasks, no failures or aborts.
- HTML reports were produced for RNA01, RNA02 and unassigned.

Fallback reads can still contain adapters and barcodes.
Read retention does not establish complete adapter removal or validate
barcode assignments.

## Pooled workflow validation

Run: reverent_wegener, completed on 2026-10-07 on MOGON.
With demultiplexing disabled, the workflow completed 18 successful tasks,
with no failures or aborted tasks. The pooled HTML report was created.
Read-retention counts matched the demultiplexed run.

Four small helper tests also passed: no trimmed reads, all reads trimmed,
mixed output, and accepted reads absent from trimmed output.
Fallback FASTQ records remained unchanged in these tests.

## Experimental fallback context trimming

The optional `barbell_trim_fallback_context` parameter defaults to false.
Enable it in the YAML parameter file:

    barbell_trim_fallback_context: true

This option applies only to fallback reads. It requires DRB004 query
sequences with the supported flanks. A read is trimmed before the left
flank only when there is exactly one exact flank-plus-18-base-barcode-prefix
candidate in the entire read, at least 10 immediately preceding A bases,
and at most 250 bases from the flank start to the read end.
Sequence and quality are cut together; the preceding poly-A sequence is
retained. These reads remain unassigned.

In the tested dataset, 235,214 of 964,873 fallback reads were trimmed;
729,659 were unchanged. All 1,107,032 input reads were retained.
The standalone experiment removed 7,807,178 bases and preserved all
51,107 fallback rRNA hits, their reference intervals and orientations.

The full pooled workflow `sick_knuth` completed with 17 successful tasks,
1 cached task, no failures or aborted tasks, and an HTML report.
Compared with the pooled run without context trimming:
- All 57,667 rRNA-aligned read IDs were identical.
- All 57,667 per-read poly-A estimates were identical.
- Modification quantification was identical at all 6,858 reference positions.

These results apply to the tested DRB004 yeast dataset. They do not prove
complete adapter removal, barcode accuracy, or equivalence to Porechop
for other datasets. The option remains experimental and disabled by default.

In retention_counts.tsv, fallback_raw counts all fallback reads, including
context-trimmed reads when enabled. fallback_context_trimmed counts those
trimmed by this option; fallback_unmodified counts those left unchanged.
