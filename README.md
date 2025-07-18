

## Project structure
```plaintext
HM_Forest/
├── DataTypes.h                      # Defines data structures and types
├── GeotiffRead.cpp                  # Handles reading GeoTIFF files
├── GeotiffWrite.cpp                 # Handles writing GeoTIFF files
├── HMFFIST.cpp                      # Main program
├── README.md                        # Project documentation
└── Data/
    └── Input/
        └── Reach_2/
            ├── 2_parameter.txt                  # Parameter configuration file
            ├── config.txt                       # Configuration file
            ├── Reach_2.zip                      # Compressed input dataset:
            │   ├── DEM_Reach_2_2m_OPTIMALfelt.tif      # Digital Elevation Model
            │   ├── River Segment Full OPTIMAL.tif      # Segmented river stream
            │   ├── src_dir_Reach_2.tif                 # Flow direction raster
            │   ├── Unet_Reach_2_v2.tif                 # Flood probability map (UNet output)
            │   ├── Reach_2.tif                         # Original image

```


## Setup & Run

**Create conda environment**
```shell
conda create -p /ENVIRONMENT/PATH/HMF
```


**Install gdal**
```shell
conda install -c conda-forge/label/cf202003 gdal
```


**Compile the code**
```shell
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/ENVIRONMENT/PATH/HMF/lib
g++ -std=c++11 -lgdal -I /ENVIRONMENT/PATH/HMF/include -L /ENVIRONMENT/PATH/HMF/lib HMFFIST.cpp GeotiffRead.cpp GeotiffWrite.cpp -O3 -o HMFFIST
```

**Run the code**
```shell
./HMFFIST ./Data/Input/Reach_2/config.txt
```


## Configuration

Create a `config.txt` file and set the following parameters:

```plaintext
./Data/Input/Reach_2/                           #Input folder
src_dir_Reach_2.tif                             #Input direction file
Unet_Reach_2_v2.tif                             #Input probability file. 
DEM_Reach_2_2m_OPTIMALfel.tif                   #Input elevation file
2_parameter.txt                                 #Input parameters file
River Segment Full OPTIMAL.tif                  #Input segmented river stream file
./Data/Result/                                  #Output folder
./Data/Result/Reach_2/                          #Output location for a test case
TC1_Prediction_HMTFIST_Unetv2.tif               #Output prediction file without regularization
TC1_Prediction_HMTFIST_Unetv2.txt               #Output as list
TC1_Prediction_HMTFIST_Unetv2_Viterbi.tif       #Output prediction file with regularization
TC1_Prediction_HMTFIST_Unetv2_Viterbi.txt       #Output as list
TC1                                             #Test case name
```
We also need to create a `2_parameter.txt` to config the following parameters:
```plaintext
2                   #Test case ID
0.001               #Epsilon
0.4                 #Pi
0.2                 #This parameters are not used in this program but required
0.5                 #This parameters are not used in this program but required
0                   #This parameters are not used in this program but required
1                   #This parameters are not used in this program but required
0.1                 #Lambda
1000                #num of intervals
```
Please remove all comments in `config.txt` and `parameters.txt`.
## Citation

If you find this code useful, please cite our paper:
```BibTex
@inproceedings{jiang2023hidden,
  title={A Hidden Markov Forest Model for Terrain-Aware Flood Inundation Mapping from Earth Imagery},
  author={Jiang, Zhe and Zhang, Yupu and Adhikari, Saugat and Yan, Da and Sainju, Arpan Man and Jia, Xiaowei and Xie, Yiqun},
  booktitle={Proceedings of the 2023 SIAM International Conference on Data Mining (SDM)},
  pages={316--324},
  year={2023},
  organization={SIAM}
}
```