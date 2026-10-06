# Dorado 2.1.2 container

Linux x86-64 image with Dorado 2.1.2, samtools, minimap2 and pod5.
The workflow default remains Dorado 0.7.2; select this image explicitly.

## Build

Requires Python 3, wget, tar and support for unprivileged Apptainer builds.
Choose a build directory outside the repository.

```bash
repo_dir=/path/to/wf-nanoribolyzer
build_dir=/path/to/container-build
mkdir -p "$build_dir"
cd "$build_dir"
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-2.1.2-linux-x64.tar.gz
tar -xzf dorado-2.1.2-linux-x64.tar.gz
python3 "$repo_dir/containers/prepare_dorado212.py" \
    "$build_dir/dorado-2.1.2-linux-x64" "$build_dir/prepared"
cd "$build_dir/prepared"
singularity build "$build_dir/dorado212_nanoribolyzer.sif" \
    "$repo_dir/containers/dorado212.def"
singularity exec --cleanenv "$build_dir/dorado212_nanoribolyzer.sif" \
    bash -e -c 'dorado --version; command -v samtools; command -v minimap2; command -v pod5'
```

Preparation fixes cuDNN links in a separate copy, preserving the source.

## Use

Add to your YAML parameters:
```yaml
dorado_container: /path/to/dorado212_nanoribolyzer.sif
dorado_executable: /opt/dorado-2.1.2-linux-x64/bin/dorado
```

Use `-profile singularity` and configure the `dorado_basecaller` label
with `containerOptions = '--nv'` for GPU access.
Bind local model directories if required; models are not included.
Remove any site-specific container override for this label when using
`dorado_container`. Keep your SLURM resource settings.

### GPU and model example

The version check above should report Dorado 2.1.2.
It checks startup and helper availability; the workflow validation below
also exercised GPU basecalling.

The following paths are examples. Replace them with your actual paths.
Example model settings matching the RNA004 models used in validation:

```yaml
dorado_models_directory: /path/to/dorado_models
dorado_rna_model: /path/to/dorado_models/rna004_sup@v6.0.0
dorado_modified_bases_models: /path/to/dorado_models/rna004_sup@v6.0.0_inosine_m6A_2OmeA@v1,/path/to/dorado_models/rna004_sup@v6.0.0_pseU_2OmeU@v1
```

Download models before running on compute nodes without internet access.
The model directory must be accessible on the compute nodes.

Add this block to your existing site configuration, replacing both
model-directory paths:

```groovy
process {
    withLabel:dorado_basecaller {
        containerOptions = '--nv --bind /path/to/dorado_models:/path/to/dorado_models'
    }
}
```

Keep your executor, queue, account and resource settings.
Remove any old Dorado host-PATH export from this label's beforeScript.
Load your site's Apptainer and Nextflow modules before launching.
Use `-profile singularity`, `-c` with your site configuration and
`-params-file` with your completed YAML.

For Barbell demultiplexing, set `preprocessing_method: barbell`,
`demultiplex: true`, and optionally select `barcodes` in the YAML.
See `references/config_template.yaml` for the remaining parameters.

## Validation

The image completed both RNA01 and RNA02 analyses on MOGON using
Dorado 2.1.2, Barbell 0.3.3 and the configured RNA004 models.
Both HTML reports were created.
