import sys,os,glob
from PIL import Image

if len(sys.argv) < 4: 
    print '**Usage: %s <folder_1> <folder_2> <new_folder> [EXT: jpg]' %(sys.argv[0])
    sys.exit(1)

folder1 = sys.argv[1]
folder2 = sys.argv[2]
foldero = sys.argv[3]
if len(sys.argv) == 5: ext = sys.argv[4]
else:                  ext = 'jpg'

if not os.path.isdir(foldero):
    os.mkdir(foldero)

sequence1 = sorted(glob.glob('%s/*.%s' %(folder1, ext)))
sequence2 = sorted(glob.glob('%s/*.%s' %(folder2, ext)))

length = min(len(sequence1), len(sequence2))

for ii in range(length): 
    print 'composite:', sequence1[ii], sequence2[ii]

    images = map(Image.open, [sequence1[ii], sequence2[ii]])
    widths, heights = zip(*(i.size for i in images))
    
    total_width = sum(widths)
    max_height = max(heights)
    
    new_im = Image.new('RGB', (total_width, max_height))
    
    x_offset = 0
    for im in images:
        new_im.paste(im, (x_offset,0))
        x_offset += im.size[0]
   
    prefix = sequence1[ii].split('/')[-1][:-4]
    new_im.save('%s/%s.%s' %(foldero, prefix, ext))
