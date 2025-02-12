from pathlib import Path

from myterial import orange
from rich import print

from brainrender import Scene

print(f"[{orange}]Running example: {Path(__file__).name}")

# Create a brainrender scene
scene = Scene(title="brain regions", atlas_name="allen_human_500um")

# Add brain regions
scene.add_brain_region("FGM")

# You can specify color, transparency...

# Render!
scene.render()
