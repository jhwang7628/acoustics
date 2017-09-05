#!/bin/bash

#./render.py

#./rename.py

if [ 0 -eq 1 ]
then
echo "padding!"
cd frames
cp frame_060.png frame_061.png
cp frame_060.png frame_062.png
cp frame_060.png frame_063.png
cp frame_060.png frame_064.png
cp frame_060.png frame_065.png
cp frame_060.png frame_066.png
cp frame_060.png frame_067.png
cp frame_060.png frame_068.png
cp frame_060.png frame_069.png
cp frame_060.png frame_070.png
cp frame_060.png frame_071.png
cp frame_060.png frame_072.png
cp frame_060.png frame_073.png
cp frame_060.png frame_074.png
cp frame_060.png frame_075.png
cp frame_060.png frame_076.png
cp frame_060.png frame_077.png
cp frame_060.png frame_078.png
cp frame_060.png frame_079.png
cp frame_060.png frame_080.png
cp frame_060.png frame_081.png
cp frame_060.png frame_082.png
cp frame_060.png frame_083.png
cp frame_060.png frame_084.png
cp frame_060.png frame_085.png
cp frame_060.png frame_086.png
cp frame_060.png frame_087.png
cp frame_060.png frame_088.png
cp frame_060.png frame_089.png
cp frame_060.png frame_090.png
cd ..
fi

ffmpeg -r 60 -i frames/frame_%03d.png -c:v libx264 -crf 0 -r 60 movie0.mp4
ffmpeg -r 60 -i frames/frame_%03d.png -c:v libx264 -crf 10 -r 60 movie1.mp4
ffmpeg -r 60 -i frames/frame_%03d.png -c:v libx264 -crf 13 -r 60 movie2.mp4
ffmpeg -r 60 -i frames/frame_%03d.png -c:v libx264 -crf 16 -r 60 movie3.mp4
ffmpeg -r 60 -i frames/frame_%03d.png -c:v libx264 -crf 19 -r 60 movie4.mp4
ffmpeg -r 60 -i frames/frame_%03d.png -c:v libx264 -crf 22 -r 60 movie5.mp4

ffmpeg -i movie0.mp4 -i ../point_1.wav -c:v copy -c:a aac -strict experimental -map 0:v:0 -map 1:a:0 output0.mp4

ffmpeg -i movie1.mp4 -i ../point_1.wav -c:v copy -c:a aac -strict experimental -map 0:v:0 -map 1:a:0 output1.mp4

ffmpeg -i movie2.mp4 -i ../point_1.wav -c:v copy -c:a aac -strict experimental -map 0:v:0 -map 1:a:0 output2.mp4

ffmpeg -i movie3.mp4 -i ../point_1.wav -c:v copy -c:a aac -strict experimental -map 0:v:0 -map 1:a:0 output3.mp4

ffmpeg -i movie4.mp4 -i ../point_1.wav -c:v copy -c:a aac -strict experimental -map 0:v:0 -map 1:a:0 output4.mp4

ffmpeg -i movie5.mp4 -i ../point_1.wav -c:v copy -c:a aac -strict experimental -map 0:v:0 -map 1:a:0 output5.mp4
