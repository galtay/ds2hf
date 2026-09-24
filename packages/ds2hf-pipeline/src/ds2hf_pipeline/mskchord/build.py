"""Build the HF repository tree from the pinned study tarball.

    <out_dir>/
      README.md                  dataset card (written by `card.write_card`)
      LICENSE                    the study's LICENSE file, verbatim
      source/msk_chord_2024.tar.gz   the tarball, byte for byte
      <config>/data.parquet      one per `data_*` file in the tarball

Shipping the tarball keeps every original byte next to its conversion: the
text formatting numbers lose, the `meta_*` files and the case lists.
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path

from ds2hf_pipeline.mskchord import convert, source


@dataclass(frozen=True)
class ConfigInfo:
    name: str
    member: str
    rows: int
    columns: int


def expected_files(configs: list[str]) -> set[Path]:
    """Every file a clean build holds; anything else would be published by accident."""
    return {
        Path("README.md"),
        Path("LICENSE"),
        Path("source") / source.TARBALL_NAME,
        *(Path(config) / "data.parquet" for config in configs),
    }


def build(tarball: Path, out_dir: Path) -> list[ConfigInfo]:
    got = source.sha256(tarball)
    if got != source.TARBALL_SHA256:
        raise ValueError(f"{tarball} sha256 is {got}, pinned {source.TARBALL_SHA256}")
    members = source.read_members(tarball)

    infos = []
    for member in sorted(members):
        config = convert.config_name(member)
        if config is None:
            continue
        table = convert.read_table(member, members[member])
        convert.write_table(table, out_dir / config / "data.parquet")
        infos.append(ConfigInfo(config, member, table.num_rows, table.num_columns))

    (out_dir / "LICENSE").write_bytes(members["LICENSE"])
    (out_dir / "source").mkdir(parents=True, exist_ok=True)
    shutil.copyfile(tarball, out_dir / "source" / source.TARBALL_NAME)
    return infos
