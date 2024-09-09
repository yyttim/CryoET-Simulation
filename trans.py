import mrcfile
import numpy as np


def convert_map_to_mrc(map_file, mrc_file):
    with open(map_file, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)

        # 根据数据大小计算可能的维度
        total_size = data.size
        possible_dims = []
        for i in range(1, int(total_size ** (1 / 3)) + 1):
            if total_size % (i * i) == 0:
                j = total_size // (i * i)
                if i * i * j == total_size:
                    possible_dims.append((i, i, j))
                    possible_dims.append((i, j, i))
                    possible_dims.append((j, i, i))

        # 打印可能的维度
        print(f"Possible dimensions: {possible_dims}")

        # 假设第一个可能的维度是正确的
        nx, ny, nz = possible_dims[0]

        data = data.reshape((nz, ny, nx))

    with mrcfile.new(mrc_file, overwrite=True) as mrc:
        mrc.set_data(data)


if __name__ == "__main__":
    map_file_path = "root/autodl-tmp/temp/human/3j0g.map"
    mrc_file_path = "root/autodl-tmp/temp/human/3j0g.mrc"

    convert_map_to_mrc(map_file_path, mrc_file_path)
    print("Conversion completed successfully!")
