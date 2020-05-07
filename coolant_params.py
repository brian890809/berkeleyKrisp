

### Library of Coolant params
import math as ma
import pandas as pd
import numpy as np
import scipy.constants
import csv

#Dowtherm A
class dowtherm_A():
	def __init__(self, T):
		self.temp = T
		self.cp = np.add(1518, np.multiply(2.82,T))
		self.miu = np.divide(0.13, np.power(T,1.072))
		self.k = np.subtract(0.142, np.multiply(0.00016, T))
		self.rho = np.subtract(1078, np.multiply(0.85, T))

#Flibe
class flibe():
	def __init__(self, T):
		self.temp = T
		self.cp = 2415.78
		self.miu = np.divide(np.multiply(4.638, np.power(10, 5)), np.power(T, 2.79))
		self.k = np.add(0.7662, np.multiply(0.0005, T))
		self.rho = np.subtract(2279.92, np.multiply(0.488, T))

#FLiNaK
class flinak():
	def __init__(self, T):
		T = np.add(T, 273.15)
		self.temp = T
		self.cp = 1905.57
		if T < float(970) and T > float(770):
			self.miu = np.multiply(np.multiply(2.487, np.divide(1, np.power(10, 4))), np.exp(np.divide(4478.62, T)))
		else:
			print("T out of bound for miu")
			assert False
		if T < float(1080) and T > float(790):
			self.k = np.add(0.36, np.multiply(5.6, np.multiply(np.divide(1, np.power(10, 4)), T)))
		else:
			print("T out of bound for k") 
			assert False
		if T < float(1170) and T > float(940):
			self.rho = np.subtract(2729.3, np.multiply(0.488, T))
		else:
			print("T out of bound for rho")
			assert False
