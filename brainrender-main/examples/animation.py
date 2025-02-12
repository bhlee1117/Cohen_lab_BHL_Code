"""
This example shows how to create an animated video by specifying
the camera parameters at a number of key frames
"""

from pathlib import Path

from myterial import orange
from rich import print

from brainrender import Animation, Scene

print(f"[{orange}]Running example: {Path(__file__).name}")


# Create a brainrender scene
scene = Scene(title="brain regions", inset=False)

# Add brain regions
scene.add_brain_region("TH")

anim = Animation(scene, Path.cwd(), "brainrender_animation")

# Specify camera position and zoom at some key frames
# each key frame defines the scene's state after n seconds have passed
anim.add_keyframe(0, camera="top", zoom=1)
anim.add_keyframe(1.5, camera="sagittal", zoom=0.95)
anim.add_keyframe(3, camera="frontal", zoom=1)
anim.add_keyframe(4, camera="frontal", zoom=1.2)

# Make videos
anim.make_video(duration=5, fps=15, resetcam=True)
