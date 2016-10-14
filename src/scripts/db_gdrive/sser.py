import sys
import wx
import time
import os
from paraview.simple import *

if __name__ == '__main__':
    model = sys.argv[1]
    curr_path = os.getcwd()

    app = wx.App()
    screen = wx.ScreenDC()

    OpenDataFile(os.path.join(curr_path, 'temp')+ '/'  + model + '.obj')
    Show()
    Render()
    # always comes up at the same place on the screen, top left corner

    time.sleep(0.5)

    bmp = wx.EmptyBitmap(400, 400)
    mem = wx.MemoryDC(bmp)
    mem.Blit(0, 0, 400, 400, screen, 2, 23)
    del mem
    # save it in imgs and move back just in case
    os.chdir(os.path.join(curr_path, 'imgs'))
    bmp.SaveFile(model + '.png', wx.BITMAP_TYPE_PNG)
    os.chdir(curr_path)

    Delete()

