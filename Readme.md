# Pop Rec 
PopRec is an image processing pipeline, tailoured for fully automated resonstructions of spare labeled neurons from volumetric datasets
This Pipeline is composed from 3 tools each with its own Matlab GUI for easy impelementation and hopefuly compatebility: 
1) Somas to Cubes - This tool takes a 3D image stack, detects Somata by intensity and size and augments the data into a 3D cubes with a defined size in pixel range, this also outputs the coordinates of the cubes and all the somata shared by that volume.  
2) Cubes to Trees - Reconstructs the center neuorn while acconting for surrounding cells in the same cube. This tool is based of the Trees toolbox by Herman Cuntz (Ref). By selecting a range of parameters the alrorithm aoutputs different versions of that neuons for later inspection and quality control.  
3) REconstruction validation 
This app is used to review the reconstructed neurons in a randomised fashion and the experiemnter should rank each cell according to the reconstruction quality, This is can then be used to train a classifier to rank each cell and predict which reconstructions are good and which are not 

