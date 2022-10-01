# Pop Rec 
PopRec is an image processing pipeline, tailored for fully automated reconstruction of spare labelled neurons from volumetric datasets
This Pipeline is composed from 3 tools each with its own Matlab GUI for easy implementation and hopefully compatibility: 
1) Somas to Cubes - This tool takes a 3D image stack, detects Somata by intensity and size and augments the data into a 3D cubes with a defined size in pixel range, this also outputs the coordinates of the cubes and all the somata shared by that volume.  
2) Cubes to Trees - Reconstructs the center neuron while accounting for surrounding cells in the same cube. This tool is based of the Trees toolbox by Herman Cuntz (Cuntz H, Forstner F, Borst A, Häusser M (2010). One rule to grow them all: A general theory of neuronal branching and its practical application. PLoS Comput Biol 6(8): e1000877.). By selecting a range of parameters the algorithm outputs different versions of that neurons for later inspection and quality control.  
3) Reconstruction validation 
This app is used to review the reconstructed neurons in a randomised fashion and the experimenter should rank each cell according to the reconstruction quality, This is can then be used to train a classifier to rank each cell and predict which reconstructions are good and which are not 

In addition to this app you can also find other functions and notebooks used for analysis in this project and also Intensify3D+ which is build to handle normalisation across biological samples 

