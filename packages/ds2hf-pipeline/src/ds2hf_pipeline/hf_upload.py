from __future__ import annotations

from pathlib import Path

from huggingface_hub import HfApi


def upload_dataset(
    processed_dir: Path,
    repo_id: str,
    private: bool = True,
    commit_message: str | None = None,
    token: str | None = None,
    delete_patterns: list[str] | None = None,
) -> str:
    """Push <processed_dir> to a HuggingFace dataset repo, creating it if needed.

    Shared by every source. The directory is a built dataset tree: a
    README.md with `configs:` YAML frontmatter and the files it declares, as
    a source's `build` command writes it.

    `delete_patterns` lets the upload remove existing files in the remote
    repo that match a glob (e.g. `["*/train.parquet"]` when migrating away
    from an old layout). Files are deleted as part of the same commit as the
    upload, so the repo never appears in a half-migrated state.

    Auth: pass `token` explicitly, or set HF_TOKEN, or run `huggingface-cli login`
    first. Returns the URL of the dataset on success.
    """
    if not processed_dir.exists():
        raise FileNotFoundError(
            f"{processed_dir} does not exist. Run the source's build command first."
        )
    if not (processed_dir / "README.md").exists():
        raise FileNotFoundError(
            f"{processed_dir}/README.md missing. Re-run the source's build command."
        )

    api = HfApi(token=token)
    api.create_repo(
        repo_id=repo_id,
        repo_type="dataset",
        private=private,
        exist_ok=True,
    )
    # `create_repo` applies `private` only when it actually creates the repo:
    # with `exist_ok=True` an existing repo makes the whole call a no-op and
    # the flag is silently ignored. That turned a `--public` push to a repo
    # first created as a private staging copy into a private one, reported as
    # success. Assert the visibility we were asked for on every push instead.
    if api.repo_info(repo_id=repo_id, repo_type="dataset").private != private:
        api.update_repo_settings(repo_id=repo_id, repo_type="dataset", private=private)
    api.upload_folder(
        folder_path=str(processed_dir),
        repo_id=repo_id,
        repo_type="dataset",
        commit_message=commit_message or "Update dataset",
        delete_patterns=delete_patterns,
    )
    return f"https://huggingface.co/datasets/{repo_id}"
