#!/usr/bin/env python

from os import listdir, rename, chdir

chdir('frames');

i = 0
for x in sorted(listdir('.')):
  if not x.endswith('.png'):
    continue
  newname = 'frame_' + str(i).zfill(3) + '.png'
  rename(x, newname)
  i += 1
