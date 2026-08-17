"""Build a measurement-only pipeline from the original CellProfiler pipeline."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


LOAD_DATA = """LoadData:[module_num:1|svn_version:'Unknown'|variable_revision_number:6|show_window:False|notes:['Load one GFAP plane and the matching Cellpose integer labels.']|batch_state:array([], dtype=uint8)|enabled:True|wants_pause:False]
    Input data file location:Default Input Folder|
    Name of the file:load_data.csv
    Load images based on this data?:Yes
    Base image location:None|
    Process just a range of rows?:No
    Rows to process:1,100000
    Group images by metadata?:No
    Select metadata tags for grouping:None
    Rescale intensities?:No
"""


def _module_blocks(pipeline_text: str) -> dict[int, str]:
    pattern = re.compile(
        r"(?ms)^([A-Za-z][A-Za-z0-9]+):\[module_num:(\d+)\|.*?"
        r"(?=^[A-Za-z][A-Za-z0-9]+:\[module_num:|\Z)"
    )
    return {int(match.group(2)): match.group(0).rstrip() for match in pattern.finditer(pipeline_text)}


def build_adapter(source: Path, destination: Path) -> None:
    """Preserve exact measurement settings while replacing segmentation."""
    text = source.read_text(encoding="utf-8-sig")
    blocks = _module_blocks(text)
    required = {12: "MeasureObjectIntensity", 13: "MeasureObjectSizeShape", 16: "ExportToSpreadsheet"}
    missing = [number for number in required if number not in blocks]
    if missing:
        raise ValueError(f"Source pipeline is missing modules: {missing}")

    measurement_blocks: list[str] = []
    for new_number, source_number in enumerate((12, 13, 16), start=2):
        block = re.sub(
            rf"module_num:{source_number}\|",
            f"module_num:{new_number}|",
            blocks[source_number],
            count=1,
        )
        block = block.replace("Astrocytes_with_nuclei", "Cellpose_masks")
        if source_number == 16:
            block = block.replace(
                "Export all measurement types?:No",
                "Export all measurement types?:Yes",
            )
            block = block.replace(
                "Overwrite existing files without warning?:No",
                "Overwrite existing files without warning?:Yes",
            )
        measurement_blocks.append(block)

    header = (
        "CellProfiler Pipeline: http://www.cellprofiler.org\n"
        "Version:5\n"
        "DateRevision:428\n"
        "GitHash:\n"
        "ModuleCount:4\n"
        "HasImagePlaneDetails:False\n\n"
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        header + "\n\n".join([LOAD_DATA.rstrip(), *measurement_blocks]) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    build_adapter(args.source, args.destination)


if __name__ == "__main__":
    main()
