# ASR： Automated segmentation and recognition of C.elegans whole-body cells

![pipeline](https://github.com/reaneyli/ASR/blob/main/pipeline.png)


We report an automatic pipeline, ASR, which automatically segments and recognizes all 558 cells of new born L1-stage larvae. For segmentation, we propose a deep learning algorithm for cell segmentation based on distance-vector. Within each distance-vector map, pixels between separate cells have a significant difference, which is very favorable for dense cell segmentation. In recognition,  we propose an iterative alignment-based assignment recognition algorithm based on the statistical atlas. Finally, ASR obtained an AP (averaged over Intersection-over-Union thresholds) of 0.8781 on 116 L1-stage larvae (AP@0.5 0.8879, AP@0.75 0.8773).

 ## Denpendencies
To run the training and inference scrips, several dependencies are required to be installed.
To facilitate the installation process,We have also prepared a docker image which contains all the required packages. If the docker environment is available on the workstation, the docker image can be accessed by running:

    docker pull fthirty/ubuntu16

 ## Run ASR
 ### Preparing the training data

To avoid losing a large number of features by excessive downsampling, we choose image blocks as input to the model, which can be done by running the following script: 

     python  process_data/get_patch.py

Generate the corresponding distance vector map by running:

    python process_data/Downsampling_transformation.py

### Training
When data is ready,we can start training ASR-segmentation network by running:

    python ARS_sementation_train.py
Parameters can be modified in ARS_sementation_train.py such as model path, data path, etc., image scaling size modification can be modified in segmentation_src/dataset.py.
the trained model will be stored at ./MODEL/ .


### Inference
The following script can be used to predict the results of new data when the model has been trained:

    python process_data/visualize_P.py

the results will be stored at ./pred/. We can merge blocks of images into one image by running:

    python merge_img.py
Finally we post-process the image and calculate the metric score by running:

    python process_data/postprocessing.py
### Recognition
To obtain the cell center of mass points directly from the cell segmentation map, you can directly run the following script to obtain the recognition results:

    ARS_recognition.m



