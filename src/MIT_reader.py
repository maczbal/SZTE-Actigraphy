import struct
import numpy as np
from datetime import datetime
import pandas as pd

def actigraph_read_binary_copy(path_to_binary_file, timezone):
    # constants to define structure of binary file (bytes)
    length_of_metadata = 5
    length_of_page = 2048
    length_of_measurement_data_on_one_page = 2028
    length_of_one_dataset = 6  # tuple of x, y, z axis data (3*2 bytes)
    length_of_one_timestamp = 6
    number_of_datasets_on_one_page = length_of_measurement_data_on_one_page // length_of_one_dataset

    # open binary file in read mode and read into array then close it
    with open(path_to_binary_file, mode='rb') as file:
        binary_data = np.fromfile(file, dtype=np.uint8)


    # get calibration coefficients
    calibration_coefficients = np.array([struct.unpack(">d", binary_data[i:i + 8])[0] for i in range(0, 48, 8)])
    # delete calibration coefficients from binary data
    binary_data = binary_data[48:]

    # get metadata array out of file array
    metadata_array = binary_data[:length_of_metadata]

    # get number of pages from metadata (1st 2 bytes of metadata), decrement is
    # needed to exclude last page (last page may not be fully written)
    number_of_memory_pages = struct.unpack(">H", metadata_array[:2])[0] - 1

    # get sampling rate from metadata (3rd byte of metadata)
    sampling_rate_lookup = {
        0x17: 1,
        0x27: 10,
        0x37: 25,
        0x47: 50,
        0x57: 100
    }
    sampling_rate = sampling_rate_lookup.get(metadata_array[2], 0)

    # get measurement interval and precision mode from metadata (4th byte of metadata)
    measurement_interval_lookup = {
        0x00: (2, False),
        0x10: (4, False),
        0x20: (8, False),
        0x30: (16, False),
        0x08: (2, True),
        0x18: (4, True),
        0x28: (8, True),
        0x38: (16, True)
    }
    measurement_interval, is_high_precision = measurement_interval_lookup.get(metadata_array[3], (0, False))

    # get id from measurement data (5th byte of metadata)
    id = metadata_array[4]

     # delete metadata from binary data, after this operation the file array
    # will contain only the pages
    binary_data = binary_data[length_of_metadata:]

    # preallocate arrays to speed up process

    binary_timestamps = np.zeros((number_of_memory_pages * length_of_one_timestamp,), dtype=np.uint8)
    binary_measurement_data_x = np.zeros((number_of_memory_pages * length_of_measurement_data_on_one_page // 3,),dtype=np.uint8)
    binary_measurement_data_y = np.zeros((number_of_memory_pages * length_of_measurement_data_on_one_page // 3,),dtype=np.uint8)
    binary_measurement_data_z = np.zeros((number_of_memory_pages * length_of_measurement_data_on_one_page // 3,),dtype=np.uint8)

    xyz_index = 0
    timestamp_index = 0
    memory_page_byte_index = 0

    # iterate through pages, while iterating collect x, y, z measurement data
    # and timestamps into specific binary arrays
    for i in range(number_of_memory_pages):
        memory_page_start_byte = memory_page_byte_index

        # collect x, y, z measurement data of actual page
        measurement_data_byte_index = memory_page_start_byte

        for j in range(number_of_datasets_on_one_page):
            measurement_data_start_byte = measurement_data_byte_index

            binary_measurement_data_x[xyz_index] = binary_data[measurement_data_start_byte]
            binary_measurement_data_x[xyz_index + 1] = binary_data[measurement_data_start_byte + 1]
            binary_measurement_data_y[xyz_index] = binary_data[measurement_data_start_byte + 2]
            binary_measurement_data_y[xyz_index + 1] = binary_data[measurement_data_start_byte + 3]
            binary_measurement_data_z[xyz_index] = binary_data[measurement_data_start_byte + 4]
            binary_measurement_data_z[xyz_index + 1] = binary_data[measurement_data_start_byte + 5]

            xyz_index += 2
            measurement_data_byte_index += 6


        binary_timestamps[timestamp_index] = binary_data[
            memory_page_start_byte + length_of_measurement_data_on_one_page]
        binary_timestamps[timestamp_index + 1] = binary_data[
            memory_page_start_byte + length_of_measurement_data_on_one_page + 1]
        binary_timestamps[timestamp_index + 2] = binary_data[
            memory_page_start_byte + length_of_measurement_data_on_one_page + 2]
        binary_timestamps[timestamp_index + 3] = binary_data[
            memory_page_start_byte + length_of_measurement_data_on_one_page + 3]
        binary_timestamps[timestamp_index + 4] = binary_data[
            memory_page_start_byte + length_of_measurement_data_on_one_page + 4]
        binary_timestamps[timestamp_index + 5] = binary_data[
            memory_page_start_byte + length_of_measurement_data_on_one_page + 5]

        timestamp_index += 6
        memory_page_byte_index += length_of_page

    # get measurement data out of binary measurement data array

    scaling_constant_for_measurement_data = 2 / (2 ** 16) * measurement_interval
    measurement_data_x = (np.frombuffer(binary_measurement_data_x, dtype=np.int16).astype(np.float64) * scaling_constant_for_measurement_data)
    measurement_data_y = (np.frombuffer(binary_measurement_data_y, dtype=np.int16).astype(np.float64) * scaling_constant_for_measurement_data)
    measurement_data_z = (np.frombuffer(binary_measurement_data_z, dtype=np.int16).astype(np.float64) * scaling_constant_for_measurement_data)
    # measurement_data = np.hstack((measurement_data_x, measurement_data_y, measurement_data_z))

    if (len(measurement_data_x)%338 == 0):
        measurement_data_x = measurement_data_x[338:-338]
    else:
        maradekx = len(measurement_data_x) % 338
        measurement_data_x = measurement_data_x[338:-maradekx]

    if (len(measurement_data_y)%338 == 0):
        measurement_data_y = measurement_data_y[338:-338]
    else:
        maradeky = len(measurement_data_y) % 338
        measurement_data_y = measurement_data_y[338:-maradeky]

    if (len(measurement_data_z)%338 == 0):
        measurement_data_z = measurement_data_z[338:-338]
    else:
        maradekz = len(measurement_data_z) % 338
        measurement_data_z = measurement_data_z[338:-maradekz]


    Sx = calibration_coefficients[0]
    Ox = calibration_coefficients[1]
    Sy = calibration_coefficients[2]
    Oy = calibration_coefficients[3]
    Sz = calibration_coefficients[4]
    Oz = calibration_coefficients[5]


    for i in range(len(measurement_data_x)):
        measurement_data_x[i] *= Sx
        measurement_data_x[i] += Ox

    for i in range(len(measurement_data_y)):
        measurement_data_y[i] *= Sy
        measurement_data_y[i] += Oy

    for i in range(len(measurement_data_z)):
        measurement_data_z[i] *= Sz
        measurement_data_z[i] += Oz

    

    timestamps = [None] * number_of_memory_pages
    timestamp_index = 0
    for i in range(number_of_memory_pages):
        t11 = struct.unpack('>H', bytes(binary_timestamps[timestamp_index:timestamp_index + 2]))[0]
        t12 = struct.unpack('>H', bytes(binary_timestamps[timestamp_index + 2:timestamp_index + 4]))[0]
        t13 = struct.unpack('>H', bytes(binary_timestamps[timestamp_index + 4:timestamp_index + 6]))[0]
        t2 = struct.unpack('>I', struct.pack('>HH', t11, t12))[0]
        t3 = struct.unpack('>Q', struct.pack('>Q', t2))[0]
        t4 = float(t13) / 65535
        t5 = float(t3) + t4
        timestamp = pd.to_datetime(t5, unit='s', origin=pd.Timestamp('1904-01-01'))
        timestamps[i] = timestamp
        timestamp_index += 6

    timestamps = pd.to_datetime(timestamps, utc=True)
    timestamps = timestamps.tz_convert('Europe/Budapest')
    formatted_timestamps = timestamps.strftime("%Y-%m-%d %H:%M:%S")


    timestamps_np = np.array([datetime.utcfromtimestamp(int(ts.timestamp())) for ts in timestamps], dtype='datetime64')
    # Az időbélyegek között eltelt idők kiszámítása másodpercekben
    first_date = (timestamps[0])
    last_date = (timestamps[-2])
    T = last_date-first_date
    time_diffs = np.diff(timestamps_np) / np.timedelta64(1, 's')
    # Az átlagos mintavételi frekvencia számítása
    avg_sampling_frequency = 1 / np.mean(time_diffs)
    real_avg_sampling_frequency = avg_sampling_frequency * (len(measurement_data_x) // len(timestamps))
    real_avg_sampling_frequency = (round(real_avg_sampling_frequency, 3))

    deltaT = (T/(len(measurement_data_x)))

    return measurement_data_x, measurement_data_y, measurement_data_z, first_date, deltaT