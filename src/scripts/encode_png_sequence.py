#!/usr/bin/python
import sys, os

if (len(sys.argv) < 6):
  print "Usage: %(a0)s <prefix> <output> <yres> <xres> <frame rate>" % \
        {"a0":sys.argv[0]};
  print "  Note: 'xres' and 'yres' should be the resolution of the input image";
  print "        These parameters do not affect the output movie resolution";
  exit();

imagePrefix = sys.argv[1];

outputMovie = sys.argv[2];

yres = int(sys.argv[3]);
xres = int(sys.argv[4]);

frameRate = int(sys.argv[5]);

bitrate = 50 * 25 * yres * xres / 256;

cmd = "mencoder \"mf://%(pattern)s*.png\" -mf fps=%(frameRate)d -o %(output)s -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=%(bitrate)d" % \
  {"pattern":imagePrefix, "frameRate":frameRate, "output":outputMovie, "bitrate":bitrate};

os.system(cmd);
