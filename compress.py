import pandas as pd
import numpy as np
from math import e
from pairing import pair, depair
from KZG10 import KZG10, Newton
from time import time
import pandas as pd
from random import randint
import bisect, operator
import json
from tqdm import tqdm
def test_proof(data: list[int], coefficients, newton: Newton, method: str, indices: list[int] = None):
    match method:
        case "raw":
            indices = [i for i in range(len(data))]
            lower_bound = indices[0]
            upper_bound = indices[-1]
            starting_range = randint(lower_bound, upper_bound - int((upper_bound - lower_bound) * 0.3) - 1)
            ending_range = starting_range + int((upper_bound - lower_bound) * 0.3)
            #print(starting_range, ending_range)
            result = []
            ind = []
            for item in indices:
                if starting_range <= item <= ending_range:
                    result.append(data[item])
                    ind.append(item)
            #print("result: ", result)
            #print("indices: ", indices)
            start = time()
            multiproof, icoeff, zpoly = newton.generate_multi_proof(coefficients, indices, result)
            return time() - start
        
        case "rle":
            lower_bound = indices[0]
            upper_bound = indices[-1]
            starting_range = randint(lower_bound, upper_bound - int((upper_bound - lower_bound) * 0.3) - 1)
            ending_range = starting_range + int((upper_bound - lower_bound) * 0.3)
            real_indices1 = []
            real_indices2 = []
            result = []
            ind = []
            temp = indices.copy()
            if starting_range not in temp:
                bisect.insort(temp, starting_range)
            if ending_range not in temp:
                bisect.insort(temp, ending_range)
            temp = temp[temp.index(starting_range)-1:temp.index(ending_range)+2]
            temp.remove(starting_range)
            temp.remove(ending_range)
            for index in temp:
                real_indices1.append(indices.index(index))
                real_indices2.append(indices.index(index))
                result.append(index)
                ind.append(data[indices.index(index)])
            # for index, value in zip(indices, data):
            #     if starting_range <= index <= ending_range:
            #         real_indices1.append(indices.index(index))
            #         real_indices2.append(data.index(value))
            #         result.append(index)
            #         ind.append(value)

            start = time()
            multiproof, icoeff, zpoly = newton.generate_multi_proof(coefficients[0], real_indices1, result)
            multiproof, icoeff, zpoly = newton.generate_multi_proof(coefficients[1], real_indices2, ind)
            return time() - start
        
        case "rle_pairing":
            indices = []
            values = []
            for value in data:
                index, value = depair(value)
                indices.append(index)
                values.append(value)
            
            lower_bound = indices[0]
            upper_bound = indices[-1]
            starting_range = randint(lower_bound, upper_bound - int((upper_bound - lower_bound) * 0.3) - 1)
            ending_range = starting_range + int((upper_bound - lower_bound) * 0.3)
            result = []
            ind = []
            temp = indices.copy()
            if starting_range not in temp:
                bisect.insort(temp, starting_range)
            if ending_range not in temp:
                bisect.insort(temp, ending_range)
            temp = temp[temp.index(starting_range)-1:temp.index(ending_range)+2]
            temp.remove(starting_range)
            temp.remove(ending_range)
            for index in temp:
                paired = pair(index, values[indices.index(index)])
                ind.append(data.index(paired))
                result.append(index)
            # for index, value in zip(indices, values):
            #     if starting_range <= index <= ending_range:
            #         paired = pair(index, value)
            #         ind.append(data.index(paired))
            #         result.append(paired)
            start = time()
            multiproof, icoeff, zpoly = newton.generate_multi_proof(coefficients, ind, result)
            return time() - start



    







def test_raw_data(data: list[float], kzg: KZG10):
    dataInt = [int(element * 10) for element in data]
    newton = kzg.choose_method(kzg.NEWTON)
    start = time()
    coefficients = newton.interpolate(data)
    ellapsed_inter = time() - start
    # print("-" * 50)
    # print("Raw data interpolation")
    # print("Data length: ", len(data))
    # print("Interpolation time (executed once): ", ellapsed_inter)
    return len(data), ellapsed_inter

def test_rle_compressed_data(data: list[int], kzg: KZG10, window_size: int, error_tolerance: float):
    indices, values, avg_err = fixed_size_window(data, window_size, error_tolerance)
    newton = kzg.choose_method(kzg.NEWTON)
    start = time()
    coefficients1 = newton.interpolate(indices)
    coefficients2 = newton.interpolate(values)
    ellapsed = time() - start
    # print("-" * 50)
    # print("Compressed data interpolation")
    # print("Data length: ", len(values))
    # print("Interpolation time (executed twice): ", ellapsed)
    return len(values), ellapsed, avg_err

def test_rle_compressed_data_with_pairing(data: list[int], kzg: KZG10, window_size: int, error_tolerance: float):
    indices, values, avg_err = fixed_size_window(data, window_size, error_tolerance)
    pairedValues = [pair(index, value) for index, value in zip(indices, values)]
    newton = kzg.choose_method(kzg.NEWTON)
    start = time()
    coefficients = newton.interpolate(pairedValues)
    ellapsed = time() - start
    # print("-" * 50)
    # print("Compressed data with pairing interpolation")
    # print("Data length: ", len(pairedValues))
    # print("Interpolation time (executed once): ", ellapsed)
    return len(pairedValues), ellapsed, avg_err
    

def fixed_size_window(data: list[int], window_size: int, error_tolerance: float):
    errors = []
    moving_averages = []
    i = 0
    df = pd.DataFrame({"data": data})
    for i in range(len(data)):

        if i < len(data) - window_size + 1:      
            window = data[i : i + window_size]
            window_average = round(sum(window) / window_size, 1)
            if i == 0:
                for j in range(window_size):
                    
                    if abs((window_average - data[j + window_size - 1]) / data[j + window_size - 1]) >= error_tolerance:
                        moving_averages.append(data[j + window_size - 1])
                        errors.append(0.0)
                    else:
                        moving_averages.append(window_average)
                        errors.append(abs((window_average - data[j + window_size - 1]) / data[j + window_size - 1]))
            else:
                if abs((window_average - data[i + window_size - 1]) / data[i + window_size - 1]) >= error_tolerance:
                    moving_averages.append(data[i + window_size - 1])
                    errors.append(0.0)
                else:
                    moving_averages.append(window_average)
                    errors.append(abs((window_average - data[i + window_size - 1]) / data[i + window_size - 1]))
        else:
            break
    df["moving_averages"] = moving_averages
    df["errors"] = errors

    indices = []
    values = []
    indices.append(0)
    values.append(moving_averages[0])
    
    for i in range(1, len(moving_averages) - 1):
        if values[-1] != moving_averages[i]:
            indices.append(i)
            values.append(moving_averages[i])

    indices.append(len(moving_averages) - 1)
    values.append(moving_averages[-1])
    values = [int(element * 10) for element in values]
    avg_err = sum(errors) / len(errors)
    return indices, values, avg_err


def entropy(labels, base=None):
    vc = pd.Series(labels).value_counts(normalize=True, sort=False)
    base = e if base is None else base
    return -(vc * np.log(vc)/np.log(base)).sum()



if __name__ == "__main__":
    kzg = KZG10()
    df = pd.read_excel("data/1DAY-Tempereature-Humidity.xlsx", sheet_name="Sheet1")
    df = df.drop(columns=['time'])
    df = abs(df)
    entropies = {}
    categories = ['ble_humi_3', 'ble_humi_4', 'ble_temp4'] # [min, mid, max]
    window_sizes = [5, 6, 7, 8, 9, 10]
    err_tols = [0.4, 0.3, 0.2]
    raw = pd.DataFrame(columns=['Category', 'Data Length', 'Interpolation Time'])
    rle = pd.DataFrame(columns=['Category', 'Window Size', 'Error Tolerance', 'Data Length', 'Interpolation Time', 'Average Error'])
    rle_pairing = pd.DataFrame(columns=['Category', 'Window Size', 'Error Tolerance', 'Data Length', 'Interpolation Time', 'Average Error'])
    methods = {'raw': raw, 'rle': rle, 'rle_pairing': rle_pairing}
    progress_bar = tqdm(total=len(categories)+ len(categories) * len(window_sizes) * len(err_tols) * 2)
    for column in df.columns:
        entr = entropy(df[column].tolist())
        entropies[column] = entr
        if column in categories:
            raw_res = test_raw_data(df[column].tolist(), kzg)
            methods['raw'].loc[len(methods['raw'])] = [column, raw_res[0], raw_res[1]]
            progress_bar.update(1)
    for window_size in window_sizes:
        for err_tol in err_tols:
            for column in df.columns:
                if column in categories:
                    rle_res = test_rle_compressed_data(df[column].tolist(), kzg, window_size, err_tol)
                    methods['rle'].loc[len(methods['rle'])] = [column, window_size, err_tol, rle_res[0], rle_res[1], rle_res[2]]
                    progress_bar.update(1)
                    if rle_res[2] == 0:
                        tqdm.write(f"category: {column} with method: rle has average error 0")
                    rle_pairing_res = test_rle_compressed_data_with_pairing(df[column].tolist(), kzg, window_size, err_tol)
                    methods['rle_pairing'].loc[len(methods['rle_pairing'])] = [column, window_size, err_tol, rle_pairing_res[0], rle_pairing_res[1], rle_pairing_res[2]]
                    progress_bar.update(1)
                    if rle_pairing_res[2] == 0:
                        tqdm.write(f"category: {column} with method: rle_pairing has average error 0")
    progress_bar.close()
    entropies = dict(sorted(entropies.items(), key=operator.itemgetter(1)))
    json.dump(entropies, open("data/entropies.json", "w"))
    with pd.ExcelWriter('data/results.xlsx') as writer:
        for method in methods:
            methods[method].to_excel(writer, sheet_name=method)
