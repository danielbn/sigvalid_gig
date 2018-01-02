# from __future__ import division

import math
import numpy as np
import scipy.fftpack
import struct

__author__ = "Eli Shafer"


def parse_data(read_data, sampling_freq = 20833.0, differentiate = 0, type = 'apus', window_shift = 700):
    if type == 'scope':
        sampling_freq = 49019.0

    data_length = len(read_data)
    if data_length % 2:
        read_data = read_data[:-1]
    data_length_int16 = data_length / 2
    data_struct = '<' + str(data_length_int16) + 'h'
    data = np.asarray(struct.unpack(data_struct, read_data))

    delta_t = 1.0 / sampling_freq
    time = np.linspace(0.0, (data_length_int16 - 1) * delta_t, data_length_int16)

    if type == 'apus':
        #Apus calculation
        acceleration_g = (data * 2.5) / 32768.0

        # Cut time to 2 to the power of something. ie cut to 1, 2, 4, 8, etc.
        max_time = np.max(time)
        power_of_two = math.log(max_time, 2)
        power_of_two = math.floor(power_of_two)
        max_time = 1 << int(power_of_two)

        time_cut = time[np.where(time < max_time)]
        N = len(time_cut)
        time = time_cut
        acceleration_g = acceleration_g[window_shift:(N + window_shift)]

    elif type == 'scope':
        #Auguscope calculation
        acceleration_g = (data / 32768.0)
        N = len(acceleration_g)

    if differentiate:
        jerk = np.gradient(acceleration_g)
    else:
        jerk = 0

    yf = scipy.fftpack.fft(acceleration_g)
    xf = np.linspace(0.0, sampling_freq / 2.0 - (1 / time[-1] + delta_t), N / 2)
    yfa = (2.0 / N * np.abs(yf[:N // 2]))

    return (time, acceleration_g, jerk), (xf, yfa)


def open_and_read_data(file_dir, file_name, sampling_freq = 20833.0, differentiate = 0, type = 'apus'):
    with open(file_dir + file_name, 'rb') as f:
        read_data = f.read()
        if read_data[:4] == '^^%\x0A':
            read_data = read_data.split("^^$$\x0A")
            read_data = read_data[1]
    return parse_data(read_data, sampling_freq, differentiate = differentiate, type = type)

file_dir='/home/irozenberg/Development/Code/SharedDevIgal/sigvalid_gig/rawdata/'
file_name='0fc03227-73b2-4e04-b2f2-67adc6b9298a'

(time, acceleration_g, jerk), (xf, yfa)=open_and_read_data(file_dir, file_name)

import matplotlib.pyplot as pyplot
pyplot.plot(time)
print 'meow'