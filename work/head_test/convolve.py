
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
		fname = "data_left_right/head_left_%05d" % i
		alldata['left'].append(readHeadFile(fname))

	for i in range(0, num):
		fname = "data_left_right/head_right_%05d" % i
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

	maxii = 0
	for pt in range(0, 100):
		left = getData('left', pt, alldata)
		right = getData('right', pt, alldata)
		if pt == 41:
			plt.plot(left)

		left_sound = numpy.convolve(sound, left, mode='same')
		right_sound = numpy.convolve(sound, right, mode='same')

		maxi = max(max([abs(w) for w in left_sound]), max([abs(w) for w in right_sound]))
		maxii = max(maxii, maxi)

	print maxii
	plotd('left', 0, alldata)

	for pt in range(0, 100):
		left = getData('left', pt, alldata)
		right = getData('right', pt, alldata)

		left_sound = numpy.convolve(sound, left, mode='same')
		right_sound = numpy.convolve(sound, right, mode='same')

		left_sound = [w/maxii for w in left_sound]
		right_sound = [w/maxii for w in right_sound]

		# plt.plot(left_sound)
		# plt.plot(right_sound)

		frames = []
		for i in range(len(left_sound)):
			bit = int(32767*right_sound[i])
			frames.append(struct.pack("<h", bit))
			bit = int(32767*left_sound[i])
			# bit = 0
			frames.append(struct.pack("<h", bit))

		dadata = ''.join(frames)

		outputFile = wave.open("pluck_stereo_%d.wav"%pt, "w")
		outputFile.setnchannels(2)
		outputFile.setsampwidth(2)
		outputFile.setframerate(44100)
		outputFile.writeframes(dadata)
		outputFile.close()

	# plt.plot([left_sound[i]-right_sound[i] for i in range(len(left))])
	# plt.plot(left_sound)
	plt.show()