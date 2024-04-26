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

def ZMPY_CLI_BatchShapeScore(PDBFileNameA, PDBFileNameB,GridWidth):
    def ZMCal(PDBFileName,GridWidth,BinomialCache, CLMCache, CLMCache3D, GCache_complex, GCache_complex_index, GCache_pqr_linear, MaxOrder, Param, ResidueBox, RotationIndex):

        [XYZ,AA_NameList]=z.get_pdb_xyz_ca(PDBFileName);
        
        [Voxel3D,Corner]=z.fill_voxel_by_weight_density(XYZ,AA_NameList,Param['residue_weight_map'],GridWidth,ResidueBox[GridWidth]);
        Dimension_BBox_scaled=Voxel3D.shape;
        
        XYZ_SampleStruct = {
            'X_sample': np.arange(Dimension_BBox_scaled[0] + 1),
            'Y_sample': np.arange(Dimension_BBox_scaled[1] + 1),
            'Z_sample': np.arange(Dimension_BBox_scaled[2] + 1)
        }
        
        [VolumeMass,Center,_]=z.calculate_bbox_moment(Voxel3D,1,XYZ_SampleStruct);
        [AverageVoxelDist2Center,MaxVoxelDist2Center]=z.calculate_molecular_radius(Voxel3D,Center,VolumeMass,Param['default_radius_multiplier'])

        ##################################################################################
        # You may add any preprocessing on the voxel before applying the Zernike moment. #
        ##################################################################################
        
        SphereXYZ_SampleStruct=z.get_bbox_moment_xyz_sample(Center,AverageVoxelDist2Center,Dimension_BBox_scaled)
        _,_,SphereBBoxMoment=z.calculate_bbox_moment(Voxel3D,MaxOrder,SphereXYZ_SampleStruct)

        [ZMoment_scaled, ZMoment_raw]=z.calculate_bbox_moment_2_zm(MaxOrder, GCache_complex, GCache_pqr_linear, GCache_complex_index, CLMCache3D, SphereBBoxMoment)

        # ZMoment_scaled[np.isnan(ZMoment_raw)]=np.nan
        
        ZM_3DZD_invariant=z.get_3dzd_121_descriptor(ZMoment_scaled)
        
        TargetOrder2NormRotate=2
        ABList_2=z.calculate_ab_rotation(ZMoment_raw, TargetOrder2NormRotate)
        ZMList_2=z.calculate_zm_by_ab_rotation(ZMoment_raw, BinomialCache, ABList_2, MaxOrder, CLMCache,s_id,n,l,m,mu,k,IsNLM_Value)
        [ZM_2,_]=z.get_mean_invariant(ZMList_2)
    
        TargetOrder2NormRotate=3
        ABList_3=z.calculate_ab_rotation(ZMoment_raw, TargetOrder2NormRotate)
        ZMList_3=z.calculate_zm_by_ab_rotation(ZMoment_raw, BinomialCache, ABList_3, MaxOrder, CLMCache,s_id,n,l,m,mu,k,IsNLM_Value)
        [ZM_3,_]=z.get_mean_invariant(ZMList_3)
        
        TargetOrder2NormRotate=4
        ABList_4=z.calculate_ab_rotation(ZMoment_raw, TargetOrder2NormRotate)
        ZMList_4=z.calculate_zm_by_ab_rotation(ZMoment_raw, BinomialCache, ABList_4, MaxOrder, CLMCache,s_id,n,l,m,mu,k,IsNLM_Value)
        [ZM_4,_]=z.get_mean_invariant(ZMList_4)
        
        TargetOrder2NormRotate=5
        ABList_5=z.calculate_ab_rotation(ZMoment_raw, TargetOrder2NormRotate)
        ZMList_5=z.calculate_zm_by_ab_rotation(ZMoment_raw, BinomialCache, ABList_5, MaxOrder, CLMCache,s_id,n,l,m,mu,k,IsNLM_Value)
        [ZM_5,_]=z.get_mean_invariant(ZMList_5)

        MomentInvariant=np.concatenate([ z[~np.isnan(z)]     for z in [ZM_3DZD_invariant,ZM_2,ZM_3,ZM_4,ZM_5]])
        
        TotalResidueWeight=z.get_total_residue_weight(AA_NameList,Param['residue_weight_map'])

        [Prctile_list,STD_XYZ_dist2center,S,K]=z.get_ca_distance_info(XYZ);

        GeoDescriptor = np.vstack((AverageVoxelDist2Center, TotalResidueWeight, Prctile_list, STD_XYZ_dist2center, S, K))
        
        return MomentInvariant, GeoDescriptor

    
    Param=z.get_global_parameter();
    
    MaxOrder=20
    
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

    MaxN=MaxOrder+1
    
    ResidueBox=z.get_residue_gaussian_density_cache(Param)

    P=z.get_descriptor_property()

    GeoScoreScaled = []  
    ZMScoreScaled  = []  
    
    for index in range(len(PDBFileNameA)):
        MomentInvariantRawA, GeoDescriptorA =ZMCal(PDBFileNameA[index],GridWidth,BinomialCache, CLMCache, CLMCache3D, GCache_complex, GCache_complex_index, GCache_pqr_linear, MaxOrder, Param, ResidueBox, RotationIndex)
        MomentInvariantRawB, GeoDescriptorB =ZMCal(PDBFileNameB[index],GridWidth,BinomialCache, CLMCache, CLMCache3D, GCache_complex, GCache_complex_index, GCache_pqr_linear, MaxOrder, Param, ResidueBox, RotationIndex)

        ZMIndex = np.concatenate((P['ZMIndex0'], P['ZMIndex1'], P['ZMIndex2'], P['ZMIndex3'], P['ZMIndex4']))
        ZMWeight = np.concatenate((P['ZMWeight0'], P['ZMWeight1'], P['ZMWeight2'], P['ZMWeight3'], P['ZMWeight4']))
        
        # Calculating ZMScore
        ZMScore = np.sum(np.abs(MomentInvariantRawA[ZMIndex] - MomentInvariantRawB[ZMIndex]) * ZMWeight)
        
        # Calculating GeoScore
        GeoScore = np.sum(P['GeoWeight'] * (2 * np.abs(GeoDescriptorA - GeoDescriptorB) / (1 + np.abs(GeoDescriptorA) + np.abs(GeoDescriptorB))))
    
        # # Calculating paper loss, not report
        # Paper_Loss = ZMScore + GeoScore
        
        # Scaled scores, fitted to shape service
        GeoScoreScaled.append((6.6 - GeoScore) / 6.6 * 100.0)
        ZMScoreScaled.append((9.0 - ZMScore) / 9.0 * 100.0)

    return GeoScoreScaled, ZMScoreScaled


def main():
    if len(sys.argv) != 3:
        print('Usage: ZMPY_CLI_BatchShapeScore PDBFileList.txt GridWidth')
        print('       This function takes a list of paired PDB structure file paths and a grid width to generate shape analysis scores.')
        print("Error: You must provide exactly one input file and a grid width.")
        sys.exit(1)

    parser = argparse.ArgumentParser(description='Process input file that contains paths to .pdb or .txt files and a grid width.')
    parser.add_argument('input_file', type=str, help='The input file that contains paths to .pdb or .txt files.')
    parser.add_argument('grid_width', type=str, help='The grid width (must be 0.25, 0.50 or 1.0)')

    args = parser.parse_args()

    input_file = args.input_file
    if not input_file.endswith('.txt'):
        parser.error("PDB file list must end with .txt")
    
    if not os.path.isfile(input_file):
        parser.error("File does not exist")

    try:
        GridWidth = float(args.grid_width)
    except ValueError:
        parser.error("GridWidth cannot be converted to a float.")
    

    if GridWidth not in [0.25, 0.50, 1.0]:
        parser.error("grid width must be either 0.25, 0.50, or 1.0")


    with open(input_file, 'r') as file:
        lines = file.readlines()

    file_list_1 = []
    file_list_2 = []
    for line in lines:
        files = line.strip().split()
        if len(files) != 2:
            print(f"Error: Each line must contain exactly two file paths, but got {len(files)}.")
            sys.exit(1)
        file1, file2 = files

        for file in [file1, file2]:
            if not (file.endswith('.pdb') or file.endswith('.txt')):
                print(f"Error: File {file} must end with .pdb or .txt.")
                sys.exit(1)
            if not os.path.isfile(file):
                print(f"Error: File {file} does not exist.")
                sys.exit(1)
        file_list_1.append(file1)
        file_list_2.append(file2)


    GeoScoreScaled, ZMScoreScaled=ZMPY_CLI_BatchShapeScore(file_list_1,file_list_2,GridWidth)

    # print('Left, the scaled score for the geometric descriptor.')
    # print('Right, the scaled score for the Zernike moments.')
    for geo_score, zm_score in zip(GeoScoreScaled, ZMScoreScaled):
        print(f'GeoScore {geo_score:.2f} TotalZMScore {zm_score:.2f}')

if __name__ == "__main__":
    main()