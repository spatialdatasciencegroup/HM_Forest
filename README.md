## Input

Prediction image(e.g., Unet prediction)

Tau-DEM direction image

Elevation image

River segmentation image

Bank of river image

## Instruction

**Create conda environment**

conda create -p PATH/HMF

**Install gdal**

conda install -c conda-forge/label/cf202003 gdal

**Compile the code**

g++ -std=c++11 -lgdal -I PATH/HMF/include -L PATH/HMF/lib HMFFIST.cpp GeotiffRead.cpp GeotiffWrite.cpp -O3 -o HMFFIST

**Run the code**

./HMFFIST ./Data/Input/Reach_2/config.txt