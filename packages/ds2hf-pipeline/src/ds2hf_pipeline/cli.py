"""`ds2hf-pipeline <source> <command>`: one sub-command group per source."""

from __future__ import annotations

import typer
from dotenv import load_dotenv

from ds2hf_pipeline.mskchord.cli import app as mskchord_app
from ds2hf_pipeline.tcga.cli import app as tcga_app

# Load .env from cwd (or any parent) on import. override=True so the project's
# .env wins over any inherited shell variable: the HF_TOKEN here is scoped to
# this project, and we don't want a stale global token to silently take over.
load_dotenv(override=True)

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="Build public cancer datasets and publish them to the HF Hub.",
)
app.add_typer(tcga_app, name="tcga")
app.add_typer(mskchord_app, name="mskchord")

if __name__ == "__main__":
    app()
