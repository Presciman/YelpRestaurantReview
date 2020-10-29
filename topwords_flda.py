#!/usr/bin/python
from wordcloud import WordCloud
import matplotlib.pyplot as plt
import sys
from operator import itemgetter
from math import exp

def main(): 
	K = int(sys.argv[1])
	Z = int(sys.argv[2])
	Y = int(sys.argv[3])
	filename = sys.argv[4]
	png_name_prefix = filename.split('.')[0]
	N = float(int(sys.argv[5])) # number of samples
	infile = file(filename+'.assign', "r")

	Z0 = Z
	for k in range(1, K): Z0 *= Y

	numtokens = 0
	numcounts = 0

	count = {}

	count = {}
	for i in range(0, Z0): count[i] = {}

	for line in infile:
		tokens = line.split()
		if len(tokens) < 3: continue

		id = tokens.pop(0)

		for token in tokens:
			parts = token.split(":")
			z = parts[-Z0:]
			word = ":".join(parts[0:len(parts)-Z0])

			numtokens += 1
			for i in range(0, Z0):
				if word not in count[i]:
					count[i][word] = 0
				count[i][word] += int(z[i]) / N
				numcounts += int(z[i])

	infile.close()

	omegaB = 0
	infile = open(filename+".omegaB", 'r')
	for line in infile:
		parts = line.strip()
		omegaB = float(parts)
	
	omegaW = {}
	infile = open(filename+".omegaW", "r")
	for line in infile:
		parts = line.strip().split()
		word = parts[0]
		w = parts[1]

		weight = float(w)
		omegaW[word] = (weight)
	infile.close()

	print ('Omega weights:')
	print ('Background (omega^0):')
	w = 0
	words = sorted(omegaW.items(), key=itemgetter(1), reverse=True)
	for word, v in words:
		print (word, v)
		w += 1
		if w >= 30: break
	print ("\n")


	omega = {}
	for k in range(0,K):
		omega[k] = {}
		infile = open(filename+".omegaZW"+str(k), "r")
		for line in infile:
			parts = line.strip().split()
			word = parts.pop(0)
	
			y = 0
			for w in parts:
				weight = float(w)
				if y not in omega[k]: omega[k][y] = {}
				omega[k][y][word] = weight
				y += 1
		infile.close()

		if k == 0: continue

		for y in range(0,Y):
			print ('Omega weights:')
			print ("Factor %d, Component %d" % (k, y))
			w = 0
			words = sorted(omega[k][y].items(), key=itemgetter(1), reverse=True)
			for word, v in words:
				print ("  "+word, v)
				w += 1
				if w >= 30: break
			print ("\n")

	prior = {}
	for x in range(Z0):
		y = {}
		for k in range(1, K):
			y[k] = int(x / pow(Y, K-k-1)) % Y

		y[0] = int(x / pow(Y, K-1))

		prior[x] = {}
		for word in omegaW:
			sum = omegaB + 1.0*omegaW[word]

			for k in range(K): sum += omega[k][y[k]][word]

			prior[x][word] = exp(sum)

			if x not in count: count[x] = {}
			if word not in count[x]: count[x][word] = 0
			# uncomment to include prior in counts
			#count[x][word] += prior[x][word]

	beta = []
	infile = open(filename+".beta", "r")
	for line in infile:
		beta.append(float(line.strip()))
	infile.close()


	for z in range(Z):
		print ('Omega weights:')
		print ("Factor 0 (Topic), Component %d\n" % (z+0))

		w = 0
		words = sorted(omega[0][z].items(), key=itemgetter(1), reverse=True)
		d = {}
		for word, v in words:
			print ("  "+word, v)
			d[word] = v
			w += 1
			if w >= 30: break

		wordcloud = WordCloud(width=900,height=500, max_words=1628,relative_scaling=1,normalize_plurals=False).generate_from_frequencies(d)
		plt.imshow(wordcloud,interpolation='bilinear')
		plt.axis("off")
		png_name = png_name_prefix + '_Component'+str(z) +'_result.png'
		plt.savefig(png_name)

		print ("\n")

		size = pow(Y,K-1)
		for xy in range(0, size):
			x = z*size + xy

			comp = [str(z)]
			for k in range(1,K):
				sizek = pow(Y,K-1-k)
				y = int(xy / sizek)
				xy = xy % sizek
				comp.append(str(y))

			print ('Sampler counts:')
			print ('Vector: %s (Beta: %f)' % (','.join(comp), beta[x]))

			w = 0
			words = sorted(count[x].items(), key=itemgetter(1), reverse=True)
			for word, v in words:
				print (word, v)
				w += 1
				if w >= 25: break
			print ("\n")

if __name__ ==  "__main__":
  main()
