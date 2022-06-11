#!/usr/bin/python

from dataclasses import dataclass
import os
from multiprocessing import Pool
import numpy as np

import communicate_with_cpp
import backend

import copy


@dataclass
class ChannelDetails:
	in_len: int
	max_out_len: int
	deletion_probability: float
	up_to: bool = False

	def output_alphabet_size(self):
		if self.up_to:
			return (2 ** (self.max_out_len+1)) - 1
		return 2 ** self.max_out_len

	def input_alphabet_size(self):
		return 2 ** self.in_len

@dataclass
class ExperimentDetails:
	experiment_path: str
	num_processors: int
	# The required accuracy of the BAA algorithm (affects the number of iterations until it is considered converged)
	accuracy: float = 0.05

	def trans_filename(self):
		return os.path.join(self.experiment_path, 'transmitted.codewords')

	def rec_filename(self):
		return os.path.join(self.experiment_path, 'received.codewords')

	def current_Q_filename(self):
		return os.path.join(self.experiment_path, 'current_Q.arr')

	def log_den_fn(self, i: int):
		return os.path.join(self.experiment_path, f'log_den_{i}.arr')

	def alpha_fn(self, i: int):
		return os.path.join(self.experiment_path, f'alpha_{i}.arr')

	def log_den_all_fn(self):
		return os.path.join(self.experiment_path, f'log_den_all.arr')



def prep_for_baa_run(cd: ChannelDetails, ed: ExperimentDetails):
	"""
	Prepares to run the requested BAA by producing the required codeword lists.
	"""
	backend.generate_codewords(False, cd.in_len, ed.trans_filename())
	backend.generate_codewords(cd.up_to, cd.max_out_len, ed.rec_filename())

def backend_compute_log_dens(params):
	return backend.compute_log_dens(*params)

def compute_log_dens(cd: ChannelDetails, ed: ExperimentDetails):
	"""
	Distributes the computation of the logs of the denominators needed for completing a step of the BAA algorithm.
	"""
	jump_size = int(np.ceil(cd.output_alphabet_size() / ed.num_processors))
	with Pool(ed.num_processors) as worker_pool:
		log_dens = np.concatenate(worker_pool.map(backend_compute_log_dens, [
									(ed.trans_filename(), ed.rec_filename(), start, start+jump_size, ed.current_Q_filename(),
										cd.deletion_probability, ed.log_den_fn(i), cd.in_len, cd.max_out_len, cd.up_to)
									for i, start in enumerate(range(0, cd.output_alphabet_size(), jump_size))
									]), axis=0)
	communicate_with_cpp.save_1d_array(log_dens, ed.log_den_all_fn())
	return log_dens


def backend_compute_alphas(params):
	return backend.compute_alphas(*params)

def compute_alphas(cd: ChannelDetails, ed: ExperimentDetails):
	"""
	Distributes the computation of the alphas, nearly completing a step of the BAA algorithm.
	"""
	jump_size = int(np.ceil(cd.input_alphabet_size() / ed.num_processors))
	with Pool(ed.num_processors) as worker_pool:
		alphas = np.concatenate(worker_pool.map(backend_compute_alphas, [
									(ed.trans_filename(), ed.rec_filename(), start, start+jump_size, ed.current_Q_filename(),
										cd.deletion_probability, ed.alpha_fn(i), cd.in_len, cd.max_out_len, cd.up_to, 
										ed.log_den_all_fn())
									for i, start in enumerate(range(0, cd.input_alphabet_size(), jump_size))
									]), axis=0)
	return alphas

def do_baa_step(initial_Q: np.ndarray, cd: ChannelDetails, ed: ExperimentDetails):
	communicate_with_cpp.save_1d_array(initial_Q, ed.current_Q_filename())
	log_dens = compute_log_dens(cd, ed)
	alphas = compute_alphas(cd, ed)
	alphas -= np.max(alphas)
	next_Q = np.exp(alphas)
	return next_Q / np.sum(next_Q)

def run_full_baa_algorithm(initial_Q: np.ndarray, cd: ChannelDetails, ed: ExperimentDetails):
	"""
	Runs the BAA algorithm, starting from some given initial distribution and continuing until the BAA bound
		shows that we are at most ed.accuracy from optimal.
	Returns the best distribution we found, how tightly the BAA bound connects it to the capacity and its rate.
	"""
	prep_for_baa_run(cd, ed)
	current_Q = copy.copy(initial_Q)
	while True:
		next_Q = do_baa_step(current_Q, cd, ed)
		distance = np.max(np.log(next_Q / current_Q))
		current_Q = next_Q
		if distance < ed.accuracy:
			return current_Q
