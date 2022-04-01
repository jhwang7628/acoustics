#!/usr/bin/python
#
# Function which returns the path to the scripts and project folders
import os

#os.path.realpath(__file__)

def scriptsPath():
    return os.path.dirname(__file__)


def projectPath():
    strPath = scriptsPath()
    return strPath[:strPath.rfind('/src/scripts')]

print(projectPath())

