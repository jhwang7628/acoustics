# This is the hackiest thing I've ever written and I'm ashamed I did,
# BUT it works!
#
# I'm so sorry.
#
# Paraview was causing a lot of issues in my old approach of opening the
# file, Show()ing it, Render()ing it, taking a screenshot of that region
# of the screen, and then Delete()ing it. For example, it would always fail
# to display any model besides the first one (no idea why). Because this
# would only happen if you were doing all this in the same script, I
# separated it into 2: one to do all the stuff above on ONE .obj, and one
# to call that script on every .obj file individually, getting around the
# weird Paraview issue.
#
# Paraview also throws some strange error:
#
#       *** Error in `python': corrupted double-linked list: (mem addr) ***
#       Aborted (core dumped)
#
# after finishing one round of the "inner" script, so just ignore that.

from subprocess import call

# change this to empty list, for now leave for testing
MODELS = ['proj.obj', 'proj\(1\).obj', 'proj\(2\).obj']
#with open('model_names.txt', 'r') as f:
#    for model in f:
#        MODELS.append(model.strip() + '.obj')

for m in MODELS:
    call('python sser.py ' + m, shell=True)

