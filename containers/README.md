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

## Validation

The image completed both RNA01 and RNA02 analyses on MOGON using
Dorado 2.1.2, Barbell 0.3.3 and the configured RNA004 models.
Both HTML reports were created.
