
import numpy
import matplotlib.pyplot as plt
import struct
import wave

def readHeadFile(filename):
	ret = {}
	ret['data'] = []
	with open(filename) as file:
		i = 0
		for line in file:
			if i == 0:
				w = line.split()
				ret['n'] = int(w[0])
				ret['time'] = float(w[1])
				i = 1
			else:
				ret['data'].append(float(line))

	return ret

def getData(side, point, alldata):
	x = [datum['data'][point] for datum in alldata[side]]
	return x

def plotd(side, point, alldata):
	time = [datum['time'] for datum in alldata[side]]
	# x = [datum['data'][point] for datum in alldata[side]]
	# plt.plot(time, x)
	print 1.0/(time[1]-time[0])

if __name__ == '__main__':

	alldata = {'left': [], 'right': []}

	num = 440

	for i in range(0, num):
		fname = "data/head_left_%05d" % i
		alldata['left'].append(readHeadFile(fname))

	for i in range(0, num):
		fname = "data/head_right_%05d" % i
		alldata['right'].append(readHeadFile(fname))

	# plt.hold(True)

	# plotd('left', 0, alldata)
	# plt.show()
	# plotd('right', 0, alldata)
	
	test = wave.open('beethoven.wav')
	frate = test.getframerate()

	sound = []

	for i in range(test.getnframes()):
		frame = test.readframes(1)

		value = struct.unpack("<h", frame)[0]
		sound.append(value)
	
	test.close()

	# plt.plot(sound)

	maxi = max([abs(w) for w in sound]) + 0.0
	sound = [w/maxi for w in sound]

	tframes = test.getnframes()
	
	print "convolve"


	print "interpolate"

	pt = 0
	i = 0

	final_left = numpy.zeros(tframes)
	final_right = numpy.zeros(tframes)
	skip = 33000
	while i < tframes:
		print i
		begin = i
		end = min(i + skip, tframes)

		next = (pt + 100 - 5)%100
		# next = pt

		left_now = getData('left', pt, alldata)
		right_now = getData('right', pt, alldata)

		left_next = getData('left', next, alldata)
		right_next = getData('right', next, alldata)

		sound_left_now = numpy.convolve(sound[begin:end], left_now, mode='same')
		sound_right_now = numpy.convolve(sound[begin:end], right_now, mode='same')

		sound_left_next = numpy.convolve(sound[begin:end], left_next, mode='same')
		sound_right_next = numpy.convolve(sound[begin:end], right_next, mode='same')

		for j in range(begin, end):
			t = (j-begin+0.0)/(end-begin + 0.0)
			interp = ((6*t - 15)*t + 10)*t*t*t
			t = interp

			final_left[j] = (1-t)*sound_left_now[j-begin] + t*sound_left_next[j-begin]
			final_right[j] = (1-t)*sound_right_now[j-begin] + t*sound_right_next[j-begin]

			# final_left[j] = sound[j]

		pt = next
		i = i + skip

	maxi = 0
	maxi = max(maxi, max([abs(w) for w in final_left]))
	maxi = max(maxi, max([abs(w) for w in final_right]))

	final_left = final_left/maxi
	final_right = final_right/maxi

	plotd('left', 0, alldata)

	print "writing"

	frames = []
	for i in range(len(final_left)):
		bit = int(32767*final_right[i])
		frames.append(struct.pack("<h", bit))
		bit = int(32767*final_left[i])
		# bit = 0
		frames.append(struct.pack("<h", bit))

	dadata = ''.join(frames)

	outputFile = wave.open("beethoven_hrtf.wav", "w")
	outputFile.setnchannels(2)
	outputFile.setsampwidth(2)
	outputFile.setframerate(44100)
	outputFile.writeframes(dadata)
	outputFile.close()

	print "finished"

	# plt.plot([left_sound[i]-right_sound[i] for i in range(len(left))])
	# plt.plot(left_sound)
	plt.show()