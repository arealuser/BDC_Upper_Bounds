#!/usr/bin/python

import subprocess
import communicate_with_cpp
import logging
BINARY_PATH = "/mnt/d/Work/Current Projects/BDC/Better Lower Bounds/BAA_in_cpp/backend/bit_channel_slave.out"
GENERATE_CODEWORDS = "gen_codewords";
COMPUTE_DENOMS = "denominators";
COMPUTE_ALPHAS = "alphas";
COMPUTE_RATE = "rate";

def run_backend(*params):
	logging.info('Backend called with the following args:'.encode('utf-8'))
	logging.info(' '.join(["'" + BINARY_PATH + "'"] + [str(x) for x in params]).encode('utf-8'))
	return subprocess.Popen([BINARY_PATH] + [str(x) for x in params])

def generate_codewords(up_to: bool, max_len: int, output_filename: str, transmitted: int):
	"""
	Uses the backend to generate a file with the codewords
	"""
	run_backend(GENERATE_CODEWORDS, int(up_to), max_len, output_filename, transmitted).wait()

def compute_log_dens(transmitted_codewords_filename: str, received_codewords_filename: str, 
	start: int, end: int, Q_array_filename: str, deletion_probability: float, output_file_name: str, 
	input_len: int, output_len: int, up_to: bool):
	"""
	Uses the backend to compute a part of the log_den arrays.
	"""
	run_backend(COMPUTE_DENOMS, transmitted_codewords_filename, received_codewords_filename, 
		start, end, Q_array_filename, deletion_probability, output_file_name, input_len, output_len, int(up_to)).wait()
	result = communicate_with_cpp.load_1d_array(output_file_name)
	return result

def compute_alphas(transmitted_codewords_filename: str, received_codewords_filename: str, 
	start: int, end: int, Q_array_filename: str, deletion_probability: float, output_file_name: str, 
	input_len: int, output_len: int, up_to: bool, log_dens: str):
	"""
	Uses the backend to compute a part of the log_den arrays.
	"""
	run_backend(COMPUTE_ALPHAS, transmitted_codewords_filename, received_codewords_filename, 
		start, end, Q_array_filename, deletion_probability, log_dens, output_file_name, input_len, output_len, int(up_to)).wait()
	result = communicate_with_cpp.load_1d_array(output_file_name)
	return result

def compute_rate(transmitted_codewords_filename: str, received_codewords_filename: str, 
	start: int, end: int, Q_array_filename: str, deletion_probability: float, output_file_name: str, 
	input_len: int, output_len: int, up_to: bool, log_dens: str):
	"""
	Uses the backend to compute a part of the log_den arrays.
	"""
	run_backend(COMPUTE_RATE, transmitted_codewords_filename, received_codewords_filename, 
		start, end, Q_array_filename, deletion_probability, log_dens, output_file_name, input_len, output_len, int(up_to)).wait()
	result = communicate_with_cpp.load_1d_array(output_file_name)
	return result[0]
