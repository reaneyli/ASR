# ASR： Automated whole-body cellular segmentation and recognition of young Caenorhabditis elegans

![pipeline](https://github.com/reaneyli/ASR/blob/main/pipeline.png)


Here, we report an automatic pipeline, ASR, which automatically segments and recognizes all 558 cells of new born L1-stage larvae. For segmentation, we propose a deep learning algorithm for cell segmentation based on distance-vector. Within each distance-vector map, pixels between separate cells have a significant difference, which is very favorable for dense cell segmentation. In recognition,  we propose an iterative alignment-based assignment recognition algorithm based on the statistical atlas. Finally, ASR obtained an AP (averaged over Intersection-over-Union thresholds) of 0.8781 on 116 L1-stage larvae (AP@0.5 0.8879, AP@0.75 0.8773).



 ## Run ASR


 ## Denpendencies
To run the training and inference scrips,several dependencies are required to be installed.
To facilitate the installation process,We have also prepared a docker image which contains all the required packages.If the docker environment is available on the workstation,the docker image can be accessed by running:
 docker pull fthirty/ubuntu16

 ### Run ASR
**Preparing the training data**

Due to the large size of the data,we need to cut the data into many pieces which can be done by running the following script: 

     python  get_patch.py

We also want to process the label so that it becomes a vector for easy training by running:

    python Downsampling_transformation.py

**Training**
When our data is ready,we can start training the 3D Unet by running:

    python train_Vector.py
Parameters can be modified in train_Vector.py such as model path, data path, etc., image scaling size modification can be modified in dataset.py.
the trained model will be stored at ./MODEL/ .


**Inference**
The test image path can be modified   in visualize_P.py.the following script is used to get the predict result:

    python visualize_P.py

the results will be stored at ./pred/. For each input image,several result files are generated.
Since the predictions are partial outcomes, we need to combine these by running:

    python merge_img.py
Finally we post-process the image and calculate the score  by running:

    python postprocessing.py





