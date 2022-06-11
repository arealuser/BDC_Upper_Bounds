#!/usr/bin/python

import numpy as np
import struct

SIZEOF_SIZE = 4
SIZEOF_FLOAT = 8

def save_array(arr: np.ndarray, filename: str):
	'''
	Given an input array and a filename, saves the array to that file in a format that the CPP code will be able to read it.
	'''
	with open(filename, 'wb') as f_out:
		f_out.write(struct.pack('II', *arr.shape))
		for row in arr:
			for val in row:
				f_out.write(struct.pack('d', val))

def load_array(filename: str):
	'''
	Given a file, loads the 2d array that is saved inside it.
	'''
	with open(filename, 'rb') as f_in:
		shape = struct.unpack('II', f_in.read(2*SIZEOF_SIZE))
		res = np.zeros(shape)
		for i in range(shape[0]):
			for j in range(shape[1]):
				res[i,j] = struct.unpack('d', f_in.read(SIZEOF_FLOAT))[0]

	return res

def save_1d_array(arr: np.ndarray, filename: str):
	'''
	Given an input 1d array and a filename, saves the array to that file in a format that the CPP code will be able to read it.
	'''
	with open(filename, 'wb') as f_out:
		f_out.write(struct.pack('I', *arr.shape))
		for val in arr:
			f_out.write(struct.pack('d', val))


def load_1d_array(filename: str):
	'''
	Given a file, loads the 1d array that is saved inside it.
	'''
	with open(filename, 'rb') as f_in:
		shape = struct.unpack('I', f_in.read(SIZEOF_SIZE))[0]
		res = np.zeros(shape)
		for i in range(shape):
			res[i] = struct.unpack('d', f_in.read(SIZEOF_FLOAT))[0]

	return res

def main():
	arr = np.array([[1., 2.], [3., 4.]])
	print(arr.shape)
	save_array(arr, "temp.array");


if __name__ == '__main__':
	main()