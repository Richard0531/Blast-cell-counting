import os
from cellpose import io, models
#import argparse
import numpy as np
import glob

#parser = argparse.ArgumentParser(description="Set the image file folder and files")
#parser.add_argument("Input", type=str)
#parser.add_argument("Output", type=str)
#args = parser.parse_args()
#input_file = args.Input
#input_dir = os.path.dirname(input_file)
#output_file = os.path.join(input_dir, args.Output)
#folder = '/home/yhsiao/pharmaco_imaging/Cellprofiler/Images/'
#file_path = '250304_121829_IMM Session/250304_133541_B1 40X PL FL Phase3/'
#file_name = 'B1_-1_3_1_Blue_003.tif'
#output_name = 'masks'
path = '/home/yhsiao/pharmaco_imaging/Cellprofiler/Images/Blast_round4_Shu/All/'
blue_files = glob.glob(os.path.join(path, '*Blue*'))
for file_path in blue_files:
    img = io.imread(file_path)
    masks, flows, styles = models.CellposeModel(model_type='cyto3').eval(img,diameter=80, channels=[0,0],normalize = {'normalize':True,'percentile':[1,99]})
    io.imsave(file_path+'flow1.tiff',flows[0][:,:,0])
    io.imsave(file_path+'flow2.tiff',flows[0][:,:,1])
    io.imsave(file_path+'flow3.tiff',flows[0][:,:,2])
#io.save_masks(img, masks, flows, output_file, tif=True,png=False)




















