# Variety_training_test_split

## **Function annotation:**

The function performs PCA decomposition and selects specific data points based on their distance from the origin point on the PCA coordinate system. This is called primary selection. Additionally, the function offers secondary data selection, which is based on the distance of every primary point to every remaining observation in the dataset.



### **Primary Selection:**

The function takes in a preprocess dataframe and performs PCA decomposition with a manually specified number of components. It then calculates the distance of each data point to the origin point on the PCA coordinate system and uses this metric to select certain data points. All observations are classified by a combination of the signs of all principal components (-1 or 1). For each class, a certain number of maximum distance points are chosen as corner points. Additionally, a certain number of minimum distance points are chosen as center points. 

Primary Selection work well for the sphere data, but poor for the skewd data, where the secondary selection come into play.

### **Secondary Selection:**

Also based on PCA space, Secondary Selection computes the distance between every possible combination of primary points and remaining observations (pairwise distance). It then applies a filter to the distances based on a user-defined threshold value. This value should be appropriate the create a speard data. Finally, gor each primary point, the observation with the minimum distance to it is selected.


### **Arguments:**

    1) df_prep: a preprocessed dataset that filters out features requiring PCA decomposition.

    2) n_components (int): number of PCs

    3) n_obs_each_corner (int): number of observation on each corner of factorial desgin.

    4) n_center_point (int): number of center points.

    5) second_select (boolean): activate secondary selection.
    
    6) distance_threshold (float): threshold used for filtering pairwise distance.
