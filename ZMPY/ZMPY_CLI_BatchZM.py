# MIT License
#
# Copyright (c) 2024 Jhih-Siang Lai
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import numpy as np
import pickle
import argparse
import os
import sys

import ZMPY as z

def ZMPY_CLI_BatchZM(PDBFileName,GridWidth=1.0, MaxOrder=6, MaxTargetOrder2NormRotate=5, Mode=0):

    Param=z.get_global_parameter()
    
    BinomialCacheFilePath = os.path.join(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache_data'), 'BinomialCache.pkl')
    with open(BinomialCacheFilePath, 'rb') as file: # Used at the entry point, it requires __file__ to identify the package location
    # with open('./cache_data/BinomialCache.pkl', 'rb') as file: # Can be used in ipynb, but not at the entry point. 
        BinomialCachePKL = pickle.load(file)

    LogCacheFilePath=os.path.join(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache_data'), 'LogG_CLMCache_MaxOrder{:02d}.pkl'.format(MaxOrder))
    with open(LogCacheFilePath, 'rb') as file: # Used at the entry point, it requires __file__ to identify the package location
    # with open('./cache_data/LogG_CLMCache_MaxOrder{:02d}.pkl'.format(MaxOrder), 'rb') as file: # Can be used in ipynb, but not at the entry point. 
        CachePKL = pickle.load(file)
    
    # Extract all cached variables from pickle. These will be converted into a tensor/cupy objects for ZMPY_CP and ZMPY_TF.
    BinomialCache=BinomialCachePKL['BinomialCache']

    # GCache, CLMCache, and all RotationIndex
    GCache_pqr_linear=CachePKL['GCache_pqr_linear']
    GCache_complex=CachePKL['GCache_complex']
    GCache_complex_index=CachePKL['GCache_complex_index']
    CLMCache3D=CachePKL['CLMCache3D']
    CLMCache=CachePKL['CLMCache']
    RotationIndex=CachePKL['RotationIndex']

    # RotationIndex is a structure, must be [0,0] to accurately obtain the s_id ... etc, within RotationIndex.
    s_id=np.squeeze(RotationIndex['s_id'][0,0])-1
    n   =np.squeeze(RotationIndex['n'][0,0])
    l   =np.squeeze(RotationIndex['l'][0,0])
    m   =np.squeeze(RotationIndex['m'][0,0])
    mu  =np.squeeze(RotationIndex['mu'][0,0])
    k   =np.squeeze(RotationIndex['k'][0,0])
    IsNLM_Value=np.squeeze(RotationIndex['IsNLM_Value'][0,0])-1

    ResidueBox=z.get_residue_gaussian_density_cache(Param)
    
    def core(PDBFileName,Mode):
        [XYZ,AA_NameList]=z.get_pdb_xyz_ca(PDBFileName)
    
        [Voxel3D,Corner]=z.fill_voxel_by_weight_density(XYZ,AA_NameList,Param['residue_weight_map'],GridWidth,ResidueBox[GridWidth])
        Dimension_BBox_scaled=Voxel3D.shape
        
        XYZ_SampleStruct = {
            'X_sample': np.arange(Dimension_BBox_scaled[0] + 1),
            'Y_sample': np.arange(Dimension_BBox_scaled[1] + 1),
            'Z_sample': np.arange(Dimension_BBox_scaled[2] + 1)
        }
        
        [VolumeMass,Center,_]=z.calculate_bbox_moment(Voxel3D,1,XYZ_SampleStruct)
    
        [AverageVoxelDist2Center,MaxVoxelDist2Center]=z.calculate_molecular_radius(Voxel3D,Center,VolumeMass,Param['default_radius_multiplier'])
    
        ##################################################################################
        # You may add any preprocessing on the voxel before applying the Zernike moment. #
        ##################################################################################
        
        SphereXYZ_SampleStruct=z.get_bbox_moment_xyz_sample(Center,AverageVoxelDist2Center,Dimension_BBox_scaled)

        _,_,SphereBBoxMoment=z.calculate_bbox_moment(Voxel3D,MaxOrder,SphereXYZ_SampleStruct)
        
        [ZMoment_scaled, ZMoment_raw]=z.calculate_bbox_moment_2_zm(MaxOrder, GCache_complex, GCache_pqr_linear, GCache_complex_index, CLMCache3D, SphereBBoxMoment)
    
        # Mode == 0 is the default, Canterakis normalisation only.
        # Mode == 1 is for 3DZD's 121 norm.
        # Mode == 2 is for both 0 and 1
        ZMList = []
        if Mode == 0:
            for TargetOrder2NormRotate in range(2, MaxTargetOrder2NormRotate+1):
                ABList = z.calculate_ab_rotation(ZMoment_raw, TargetOrder2NormRotate)
                ZM = z.calculate_zm_by_ab_rotation(ZMoment_raw, BinomialCache, ABList, MaxOrder, CLMCache,s_id,n,l,m,mu,k,IsNLM_Value)
                ZM_mean, _ = z.get_mean_invariant(ZM)
                ZMList.append(ZM_mean)
        elif Mode == 1:
            ZM_3DZD_invariant=z.get_3dzd_121_descriptor(ZMoment_scaled)
            ZMList.append(ZM_3DZD_invariant)
        elif Mode == 2:
            ZM_3DZD_invariant=z.get_3dzd_121_descriptor(ZMoment_scaled)
            ZMList.append(ZM_3DZD_invariant)
            for TargetOrder2NormRotate in range(2, MaxTargetOrder2NormRotate+1):
                ABList = z.calculate_ab_rotation(ZMoment_raw, TargetOrder2NormRotate)
                ZM = z.calculate_zm_by_ab_rotation(ZMoment_raw, BinomialCache, ABList, MaxOrder, CLMCache,s_id,n,l,m,mu,k,IsNLM_Value)
                ZM_mean, _ = z.get_mean_invariant(ZM)
                ZMList.append(ZM_mean)
        return np.concatenate([z[~np.isnan(z)] for z in ZMList])


    ZM  = []  # List to store the calculated ZM
    for f in PDBFileName:
        ZM.append(core(f,Mode))

    return ZM

def main():
    if len(sys.argv) != 6:
        print('Usage: ZMPY_CLI_BatchZM PDBFileList GridWidth MaximumOrder NormOrder Mode')
        print('    This function computes the Zernike moment based on the specified maximum order, normalization order, and voxel gridding width.')
        print("Error: You must provide exactly five input arguments.")
        sys.exit(1)

    parser = argparse.ArgumentParser(description="Process a .txt file that contains paths to .pdb or .txt files.")
    parser.add_argument('input_file', type=str, help='The input file to process (must end with .txt) containing paths to .pdb or .txt files.')
    parser.add_argument('GW', type=float, choices=[0.25, 0.50, 1.00], help='Grid width must be 0.25, 0.50 or 1.00.')
    parser.add_argument('MaxOrder', type=int, choices=[6, 20, 40], help='Maximum order of calculating ZM must be 6, 20, or 40.')
    parser.add_argument('MaxN', type=int, help='Maximum normalisation order must be an integer number.')
    parser.add_argument('Mode', type=int, choices=[0, 1, 2], help='Mode must be 0, 1 or 2.')

    args = parser.parse_args()

    input_file = args.input_file
    if not input_file.endswith('.txt'):
        parser.error("File must end with .txt")
    
    if not os.path.isfile(input_file):
        parser.error("File does not exist")

    if args.MaxN > args.MaxOrder or args.MaxN < 2:
        parser.error("Maximum normalisation order must be larger than 2 and less than or equal to the maximum order of calculating ZM.")
        

    with open(input_file, 'r') as file:
        lines = file.readlines()

    pdb_file_list=[]
    for line in lines:
        pdb_file = line.strip()

        if not (pdb_file.endswith('.pdb') or pdb_file.endswith('.txt')):
            parser.error("File must end with .pdb or .txt")
        
        if not os.path.isfile(pdb_file):
            parser.error("File does not exist")

        pdb_file_list.append(pdb_file)

    Result=ZMPY_CLI_BatchZM(pdb_file_list,args.GW,args.MaxOrder,args.MaxN,args.Mode)

    np.set_printoptions(threshold=Result[0].size)

    for x in Result:
        print(x)
        print('\n')    

if __name__ == "__main__":
    main()