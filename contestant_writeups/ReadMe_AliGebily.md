1- Introduction
----------------------------
----------------------------
The solution mainly handles the problem as a supervised machine learning problem, by predicting gene expression value using fixed number of statistical qunatites about observations. So, we needed to pre-process our data to be in tabular format. In simple steps, we do the following
1. Load input data of all measurements for all experiments, and remove outliers
2. For each measurement's observations, We create set of features by calculating statistical quantities for these observations, like count, mean, variance, min, max, median, and so on. Actually we calculate more than 100 statistical quantity representing each measurement's observations. We do this operation before training phase and prediction phase, but for training phase we save results to cvs files to save time during model training and tunning.
3. During training phase, we load all saved data, i.e statistical quantities saved in csv format, for all measurements observations from csv files, then we feed them into two xgboost models, one for trainging model to predict higher proportion gene value, and another one for predicting lower proportion gene value. Then, we save the two models to disk.
4. In prediction phase, we repeat steps 1 & 2 but without saving the calculated statistical quantities about observations to csv files. Then we load the two saved xgboost models to predict higher proportion gene and lower proportion gene values. 
5. Finally, We reformat predictions to be saved into gct format.

2- Details
----------------------------
----------------------------
Before discussing our approach for solving the problem, we will define our inputs and outputs 
- inputs
    1. observations(numbers) for the measurement
    2. the ratio of weights(proportions) of the two measured genes is `2:1`
    3. By analyzing dpk and litmus files away from each other, we discovered that data in both categories comes from two different sources, or from the same source but with different conditions or configurations. So, we consider source(dpk or litmus) while training model, and also at prediction. Actually we detect it from plate name.
- outputs: we need to identity the value of higher proportion gene(`y1`), and the value of lower propotion gene(`y0`).

Approach sequence
--------------------
- The first question for me was: is it supervised or un-supervised learning machine problem? 
Based on benchmark solution, this problem was handled as clustering(un-supervised) problem using K-means. But When we tried to run the benchmark and found how much time it takes to run, we thought that if it's handled as a supervised problem, it will give a better total score even the accuracy may decrease, because the prediction time for supervised models is lower if compared to clustering problems. So, first we decided to handle it as supervised problem, 
- So the next question came here: What are input features(columns) for each measurement?
The main issue that measurement's observations counts are different, so If we managed to find fixed number of statistical quantities that describe each measurement observations, then we could consider those quantities as our input features(columns). 

- So, we moved to next question: what are statistical qunatities that we can get from observations to keep information of observations as much a possible?

    There is a function in pandas that is usually used while data exploration called DataFrame.describe, https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.describe.html
    It gets some statistics about data samples. The pandas describe() function  returns count, mean, std and other 5 points called percentiles. Those 5 percentiles are calculated by sorting data, then return minimum, 25%(value greater than 25% of data values), 50%(median: value greater than 50% of data values),  75%(value greater than 75% of data values), and maximum.
    So we got the required statistics from observations and saved them to csv file. Then we handled the problem as a usual supervised learning problem, and built two models, one for predicting lower propotion gene value and the other for the higher proportion gene value.
    After building the model, we found that it got a reasonable accuracy and ~3 times faster than benchmark solution, so we focused on improving accuracy.
    we extended the model by including 101 percentiles from observations, not just 5 percentiles like dataframe.describe(), so we had about 104 input features(columns) including count, mean, variance, and 101 percentiles. Then we got a better score that exceeded benchmark solution, about '30.23'. For calculating percentiles, we used numpy.percentile function: https://docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.percentile.html
- We defined observations that are less than zero or observations that are too large with respect to other values to be outliers. So, we removed them.
- When we trained model, we found that it works well for large values, but its performance is poor with small values, so we decided to train model against natural logarithm of target, which finally improved performance of models, and made it consisten across all ranges.
- The next step was to try get benefit from other inputs(#2, #3) in inputs sections
    1. With respect to input #3, we added a flag called 'dpk_litmus' to identify if csv measurement record comes from dpk or litmus, and found that score improved to '31.29', as it helped our model to differentiate between measurements comes from dpk away from litmus. dkp data was more noisy than litmus data.
    2. With respect to input #2, (2:1 ratio), 
    To make better explaination, let's assume we have two separate sets of numbers with ratio `2:1`
        - `s1=[3,4,5]` => then `mean1=4`
        - `s2=[7,8,9,10,11,12]`, then `mean2=9.5`
        - `s3=[15,16,17]` => then `mean3=14` 

        If we union the sets s1 with s2 and s2 with s3, then
        - `s12=[3,4,5,7,8,9,10,11,12]`. then `mean12=7.1, median12=8`
        - `s23=[7,8,9,10,11,12,15,16,17]`. then `mean23=11.6, median23=11`

        What we can infer from those results is when the mean of lower proportion set, mean1, is less than the mean of higher proportion set, mean2, then aggreage set mean is greater than median. i.e  `mean12>median12`
        and when the mean of lower proportion set, mean3, is greater than the mean of higher proportion set, mean2, then aggreage set mean is less than median. i.e  `mean23 < median23`
        We tested this feature against our real data, and found that it's fits 90% of data, so, We added a flag column for difference between mean and median of observations that is acting as indicator if `y1 > y0`, which improved accuracy of results.

Performance notes
-------------------
- We used multicore programming to make preprocessing, training, and predictions phases run faster.
- We used `__slots__` in python instead of objects to make accessing attributes faster. You can check this: https://stackoverflow.com/questions/472000/usage-of-slots
- We scaled observation values by 100, then rounding them to integer values so that we can get benefit from faster integer calculations, and keeping data two decimal digits accurate. 

Steps to train and predict
------------------------------
- open cmd/shell and navigate to path containing source code
- build docker image
    > `docker build -t cmap/python-model .`
- run docker image. 
You need to change mounted volumes paths to match files on your machine. We mount four paths
    - path containing input data for training or prediction
    - path containing ground truth data
    - path where generated files(csvs and gct) will be saved
    - path where generated models are kept. on linux, replace `%CD%\resources:/python-model/resources` with `$(pwd)/resources:/python-model/resources`

    > `docker run --rm -it --entrypoint bash  -v D:\git\CMap-DPeak-Challenge\competitor_pack_v2\input:/input    -v D:\git\CMap-DPeak-Challenge\competitor_pack_v2\ground-truth:/ground-truth  -v D:\git\CMap-DPeak-Challenge\competitor_pack_v2\output:/output -v %CD%\resources:/python-model/resources cmap/python-model`
- generate csvs files from input txt experiments files. 
You need to provide input path containg input files, ground truth data path, path to save generated csvs files, and comma separated names of plates to generate csv for them. For each plate, we generate single csv file.

    > `python /python-model/genes-merge-percentiles-to-csv.py --inputpath /input --groundpath /ground-truth  --csvspath /output --plates DPK.CP001_A549_24H_X1_B42,LITMUS.KD017_A549_96H_X1_B42`

- train xgboost model using csv files. 
You need to specify csvspath contains csv files, and comma separated names of plates to load csv for them.
    > `python /python-model/genes-training-percentiles.py --csvspath /output --plates DPK.CP001_A549_24H_X1_B42,LITMUS.KD017_A549_96H_X1_B42`

    ### Note: generated models are saved in path relative to current running script path, in `resources` directory, and it will be loaded in prediction phase from this directory. So, you need to copy this directory as well *.py files to production server when you deploy.
- generate predictions. 
You need to specify input data path, '.gct' output file path, and name of plate to predict values for it. 
    > `python /python-model/genes-prediction.final.py --dspath /input/DPK.CP001_A549_24H_X1_B42  --out /output --create_subdir 0  --plate DPK.CP001_A549_24H_X1_B42`
- evaluate spearman correlation for generated .gct file. (optional)

    You need to specify paths of predictions, ground truth, and plate name.
    > `python /python-model/UNI-DUO-spearman.py --predictionspath /output  --groundpath  /ground-truth --plate DPK.CP001_A549_24H_X1_B42`

### Note: On machine of 2 cores & 6G RAM, all previous steps require about 70 minutes to run.

3- Questions & Answers
----------------------------
----------------------------
- Country of residence? 
    - Egypt
- Highest academic degree achieved?
    - Bachelor in Electical Engineering 
- If college educated, what was your major?
     - Communications and Electronics
- Main motivation for competing in this challenge?
    - Contributing in a project that helps getting the world better and healthier
    - Getting more experience in data science as I am new to this field.
    - Getting a prize
- How difficult/enjoyable the problem and actual competition of this MM was relative to standard MMs on Topcoder?
    - For me the problem was very enjoyable as it needed much more feature engineering to get better results. Also I learnt many machine learning topices like Gaussian Mixture models for clustering, xgboost, and python multithreading code, and so on.
    - With respect to the competition itself, I think It was a tough competition with about 294 registrants, 45 submitters, an 800 submissions. For me, It's the third MM competition to take part in effectively, but It's the first time to get a winning place, so this MM will be memorable. 

4- Feedback
----------------------------
----------------------------
- Positives:
    - Detailed explaination of challenge problem and requirements was very good point.
    - Answers to forum questions were clear
    - Having time component in score metric is very good, as inforces competitors to build a good model, not just stacking models.
    - Having automated online tester for code, not just for predictions, is very good step in topcoder
    - The distribution of prizes across 9 places was very good, as it made more users willing to compete actively which raised the bar for all competitors.
- Comments:
    - Publishing automated scorer online is delayed too much after launching challenge, so we think it's better to develop the scorer before launching the challenge. I know that challenge deadline is extended, but this increased the effort exerted during the match.
    - Because the score metric was composite including spearman correlation, AUC, and time, we were in an urgent need to get a detailed score online, so that we know which score component causes our score to go down.
    - Duration between submissions for each competitor should be defined, so that we avoid getting the online scorer down in last days.
    - Waiting many days for final announcements was not good



