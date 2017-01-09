# draw heatmaps

# if "ImportError: cannot import name imsave" is raised
# use the following command to install PIL
# sudo pip install Pillow==2.6.0

from scipy.misc import imsave
import numpy as np
import sys

if len(sys.argv) != 2:
	raise SyntaxError('Input file is missing!')

image_array = np.empty([1,120]) # 120 is the total bin number in the previous code
f_name = sys.argv[1]

f = open(f_name)
for line in f:
	line.rstrip()
	line = line.split()
	nums = [float(line[i]) for i in range(1, len(line))]
	image_array = np.append(image_array, [nums], axis=0)

imsave('heapmap.png', image_array)
