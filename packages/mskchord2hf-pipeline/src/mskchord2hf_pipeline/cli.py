from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Annotated

import typer
from dotenv import load_dotenv

from mskchord2hf_pipeline import build, card, source, verify

# override=True so the project's .env wins over any inherited shell HF_TOKEN.
load_dotenv(override=True)

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="Mirror the MSK-CHORD cBioPortal study onto the HF Hub, format conversion only.",
)

DEFAULT_DATA_DIR = Path.home() / "data" / "mskchord2hf"

DataDirOpt = Annotated[
    Path | None,
    typer.Option(
        "--data-dir",
        help="Root data dir. Defaults to $MSKCHORD2HF_DATA_DIR or $HOME/data/mskchord2hf.",
    ),
]


def _data_dir(data_dir: Path | None) -> Path:
    if data_dir is not None:
        return data_dir
    env = os.environ.get("MSKCHORD2HF_DATA_DIR")
    return Path(env).expanduser() if env else DEFAULT_DATA_DIR


def _report(checks: list[verify.Check]) -> list[str]:
    for check in checks:
        typer.echo(f"[{'PASS' if check.passed else 'FAIL'}] {check.name}: {check.summary}")
        if not check.passed:
            for line in check.details:
                typer.echo(line)
    return [c.name for c in checks if not c.passed]


@app.command("fetch")
def fetch_cmd(data_dir: DataDirOpt = None) -> None:
    """Download the study tarball into <data-dir>/raw/. sha256-pinned; safe to re-run."""
    path = source.fetch(_data_dir(data_dir) / "raw")
    typer.echo(f"{path} (sha256 {source.TARBALL_SHA256})")


@app.command("build")
def build_cmd(data_dir: DataDirOpt = None) -> None:
    """Write <data-dir>/processed/: one parquet per study data file, the tarball, the card.

    Rebuilt whole each time, so nothing stale survives.
    """
    root = _data_dir(data_dir)
    tarball = source.fetch(root / "raw")
    out_dir = root / "processed"
    if out_dir.exists():
        shutil.rmtree(out_dir)
    infos = build.build(tarball, out_dir)
    for info in infos:
        typer.echo(f"  {info.name:<30}{info.rows:>10,} rows {info.columns:>7,} cols")
    typer.echo(f"\nwrote dataset card -> {card.write_card(out_dir, infos)}")
    typer.echo("verify with: mskchord2hf-pipeline verify")


@app.command("verify")
def verify_cmd(
    data_dir: DataDirOpt = None,
    network: Annotated[
        bool, typer.Option("--network/--offline", help="Also check cBioPortal and datahub.")
    ] = True,
) -> None:
    """Check <data-dir>/processed/ cell by cell against the tarball, and the tarball upstream."""
    failed = _report(verify.verify(_data_dir(data_dir) / "processed", network=network))
    if failed:
        typer.echo(f"\n{len(failed)} check(s) failed: {', '.join(failed)}")
        raise typer.Exit(code=1)
    typer.echo("\nall checks passed")


@app.command("upload")
def upload_cmd(
    repo_id: Annotated[str, typer.Option("--repo-id", help="HF dataset repo id.")] = card.REPO_ID,
    private: Annotated[
        bool, typer.Option("--private/--public", help="Upload as private, or public.")
    ] = True,
    commit_message: Annotated[
        str | None, typer.Option("--commit-message", "-m", help="Commit message.")
    ] = None,
    data_dir: DataDirOpt = None,
) -> None:
    """Verify <data-dir>/processed/, then push it to the HF Hub.

    Private by default. Every check must pass first, including `tree_clean`,
    because `upload_folder` publishes whatever sits in the tree.
    """
    from huggingface_hub import HfApi

    processed = _data_dir(data_dir) / "processed"
    typer.echo(f"processed dir: {processed}")
    typer.echo(f"repo_id:       {repo_id} ({'private' if private else 'PUBLIC'})\n")
    failed = _report(verify.verify(processed))
    if failed:
        raise typer.BadParameter(f"not uploading: {', '.join(failed)} failed.")

    api = HfApi()
    api.create_repo(repo_id=repo_id, repo_type="dataset", private=private, exist_ok=True)
    # create_repo ignores `private` for an existing repo; assert it every push.
    if api.repo_info(repo_id=repo_id, repo_type="dataset").private != private:
        api.update_repo_settings(repo_id=repo_id, repo_type="dataset", private=private)
    api.upload_folder(
        folder_path=str(processed),
        repo_id=repo_id,
        repo_type="dataset",
        commit_message=commit_message or "Update MSK-CHORD mirror",
    )
    typer.echo(f"\nuploaded -> https://huggingface.co/datasets/{repo_id}")


if __name__ == "__main__":
    app()
