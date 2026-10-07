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
