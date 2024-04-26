# ZMPY

ZMPY: accelerating protein structure volume analysis through vectorized 3D Zernike Moments and Python-based GPU Integration

This repository is a NumPy implementation without GPU support.
For GPU support with TensorFlow and CuPy, please refer to the other two repositories:

`ZMPY_TF` supports `Tensorflow`
(https://github.com/tawssie/ZMPY_TF)

`ZMPY_CP` supports `CuPy`
(https://github.com/tawssie/ZMPY_CP)

Here presents a Python-based software package, ZMPY, to accelerate the moments computation by vectorizing the mathematical formulae, enabling their computation in graphical processing units (GPUs). The package offers popular GPU-supported libraries such as CuPy and TensorFlow along with NumPy implementations, aiming to improve computational efficiency, adaptability, and flexibility in future algorithmic development. 

## Installation

**Prerequisites:**
* ZMPY   : Python >=3.9.16, NumPy >=1.23.5
* ZMPY_CP: Python >=3.9.16, NumPy, CuPy >=12.2.0
* ZMPY_TF: Python >=3.9.16, NumPy >=1.23.5, Tensorflow >=2.12.0, Tensorflow-Probability >=0.20.1

1. Open the terminal
2. Using pip to install the package through PyPI
3. Run `pip install ZMPY` for the installation


## Usage (Google Colab notebooks)

* 3D Zernike moments: [access Colab notebook demo here]
* Shape similarity: [access Colab notebook demo here]
* Structure superposition: [access Colab notebook demo here]

## Performances

A voxel cube with dimensions of 100x100x100 was applied to perform 10,000 3D Zernike moment calculations, using 2 different maximum orders 20 and 40.
Execution times for different hardware configurations using TensorFlow, CuPy, and NumPy libraries:

### NumPy

| Order | CPU1       | CPU2       |
|-------|------------|------------|
| 20    | 33m20s     | 14m1s      |
| 40    | 951m40s    | 338m20s    |


### TensorFlow

| Order |            T4 |            RX3070Ti |            V100 |            L4 | 
|-------|---------------|---------------------|-----------------|---------------|
| 20    | 1m1s          | 0m36s               | 0m31s           | 0m39s         | 
| 40    | 24m40s        | 9m3s                | 10m54s          | 11m13s        | 

### CuPy
| Order |      T4 |      RX3070Ti |      V100 |      L4 |
|-------|---------|---------------|-----------|---------|
| 20    | 4m45s   | 2m30s         | 1m42s     | 2m50s   |
| 40    | 35m20s  | 19m19s        | 14m45s    | 18m40s  |

Note: m = minutes, s = seconds.

## Cache data for order 40

Due to GitHub's file size limitations, follow these steps to download the cache data for order 40 (1.3G) in the ZMPY package:

### 1. Locate Package Folder

- Open your terminal and execute the following command to find the folder of the ZMPY package:
- `python -c "import ZMPY; print(ZMPY.__file__)"`
- Note the path, which ends with `/User/path/ptyhon/site-packages/ZMPY/__init__.py`.

### 2. Navigate to Cache Data Folder
- Go to the `cache_data` folder at the same level as `__init__.py` file, i.e., `/User/path/ptyhon/site-packages/ZMPY/cache_data`.

### 3. Download the Cache File:
- Download the 1.3 GB max order 40 `.pkl` file to the `cache_data` folder from the link below. https://drive.google.com/uc?id=1RR1rF_5YJqaxNC5AK0Ie_8MswGb0Tttw


## Contributing

Feel free to submit pull requests for improvements or bug fixes.

************************* 


## Citation

Lai, J. S., Burley, S. K., & Duarte, J. M. (2024). ZMPY: Accelerating protein structure volume analysis through vectorized 3D Zernike moments and Python-based GPU integration. (Submitted)

## License

This project is licensed under the GNU General Public License v3.0. You can view the full license [here](https://www.gnu.org/licenses/gpl-3.0.en.html).

