#import time
import wx
from paraview.simple import *
import sys

MODELS = ['proj(1).obj', 'proj.obj']

if __name__ == '__main__':
    model = sys.argv[1]
    app = wx.App()
    screen = wx.ScreenDC()
#    size = screen.GetSize()
    for m in MODELS:
#        app = wx.App()
#        screen = wx.ScreenDC()
#        size = screen.GetSize()

        OpenDataFile('/home/vinny/Downloads/' + m)
        Show()
        Render()

        time.sleep(1)
        #        app = wx.App()

        #        screen = wx.ScreenDC()
        #        size = screen.GetSize()
        bmp = wx.EmptyBitmap(842, 552)
        mem = wx.MemoryDC(bmp)
        mem.Blit(0, 0, 842, 552, screen, 306, 206)
        del mem
        bmp.SaveFile(m[:-4] + '.png', wx.BITMAP_TYPE_PNG)

        Delete()

