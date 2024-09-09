import time
import cupy as cp
import numpy as np
from polnet import generate_synthetic_data  # 检查存在一个数据生成函数
from polnet.utils import add_noise_to_model, generate_random_shape, combine_with_random_shape

def evaluate_performance(data_func, *args, **kwargs):
    start_time = time.time()
    result = data_func(*args, **kwargs)
    cp.cuda.Stream.null.synchronize()
    end_time = time.time()
    return result, end_time - start_time

# 测试未修改前的性能
data = np.random.normal(size=(1000, 1000, 1000))
result, elapsed_time_before = evaluate_performance(np.mean, data)
print(f'Elapsed Time Before: {elapsed_time_before}s')

# 测试修改后的性能
data = cp.random.normal(size=(1000, 1000, 1000))
result, elapsed_time_after = evaluate_performance(cp.mean, data)
print(f'Elapsed Time After: {elapsed_time_after}s')

# 生成噪声和随机形状
base_model = cp.random.normal(size=(100, 3))
random_shape_verts, random_shape_faces = generate_random_shape()
enhanced_model = combine_with_random_shape(base_model, random_shape_verts)
noisy_model = add_noise_to_model(enhanced_model)
