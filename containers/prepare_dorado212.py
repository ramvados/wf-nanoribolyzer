import argparse
import shutil
from pathlib import Path

parser = argparse.ArgumentParser(description="Prepare a Dorado 2.1.2 container copy.")
parser.add_argument("source", type=Path, help="Extracted dorado-2.1.2-linux-x64 directory")
parser.add_argument("build_directory", type=Path, help="Directory for container build inputs")
args = parser.parse_args()

source = args.source.resolve()
target = args.build_directory.resolve() / "dorado-2.1.2-linux-x64"

if not (source / "bin/dorado").is_file():
    raise SystemExit("Source does not contain bin/dorado.")
if target.exists():
    raise SystemExit(f"Target already exists: {target}")

target.parent.mkdir(parents=True, exist_ok=True)
shutil.copytree(source, target, symlinks=True)

for link in sorted((target / "lib").glob("libcudnn*.so")):
    versioned = link.with_name(link.name + ".9")
    if not link.is_symlink() or not versioned.exists():
        raise SystemExit(f"Unexpected library layout: {link}")
    link.unlink()
    link.symlink_to(versioned.name)

broken = [str(p) for p in target.rglob("*") if p.is_symlink() and not p.exists()]
if broken:
    raise SystemExit("Unresolved symbolic links:\n" + "\n".join(broken))

print(f"Prepared: {target}")
