from pathlib import Path

from myterial import orange
from rich import print

from brainrender import Scene

print(f"[{orange}]Running example: {Path(__file__).name}")


# Create a brainrender scene
scene = Scene(title="slice")

# Add brain regions
th = scene.add_brain_region("TH")

# You can specify color, transparency...
mos, ca1 = scene.add_brain_region("MOs", "CA1", alpha=0.2, color="green")

# Slice actors with frontal plane
scene.slice("frontal", actors=[th])

# Slice with a custom plane
plane = scene.atlas.get_plane(pos=mos.center_of_mass(), norm=(1, 1, 0))
scene.slice(plane, actors=[mos, ca1])

# Render!
scene.render()
