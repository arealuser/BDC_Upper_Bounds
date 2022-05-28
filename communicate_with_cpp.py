#!/usr/bin/python

import numpy as np
import struct

def save_array(arr: np.ndarray, filename: str):
	'''
	Given an input array and a filename, saves the array to that file in a format that the CPP code will be able to read it.
	'''
	with open(filename, 'wb') as f_out:
		f_out.write(struct.pack('II', *arr.shape))
		for row in arr:
			for val in row:
				f_out.write(struct.pack('d', val))


def main():
	arr = np.array([[1., 2.], [3., 4.]])
	print(arr.shape)
	save_array(arr, "temp.array");


if __name__ == '__main__':
	main()