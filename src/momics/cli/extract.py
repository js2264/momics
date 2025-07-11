from pathlib import Path
import click
import cloup

from ..momics import Momics
from ..logging import logger
from .cli import cli
from .cli import Sections


@cli.command(section=Sections.io)
@click.pass_context
@click.option(
    "--tracks",
    "-t",
    type=str,
    help="Comma-separated list of track labels to restore. If not provided, all tracks will be restored.",
    default="",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(exists=False),
    required=False,
    default=".",
    help="Output directory to restore the bigwig files to (default: current directory).",
)
@cloup.argument("path", help="Path to a momics repository", metavar="MOMICS_REPO", required=True)
def extract(ctx, tracks, output, path):
    """Bulk extract all bigwig files ingested into a Momics repository."""
    m = Momics(path)

    # Select tracks to export
    tr = m.tracks()
    if tracks != "":
        tracks = tracks.split(",")
        tr = tr[tr["label"].isin(tracks)]

    # Ensure output directory exists
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)

    # Export tracks
    logger.info(f"Exporting {len(tr)} tracks to {output}")
    for _, track in tr.iterrows():
        out_f = str(Path(output, track["label"] + ".bw"))
        m.export_track(track["label"], out_f)
