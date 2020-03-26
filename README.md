# Developing a Customized Random Forest Machine Learning Approach to Identify Bacterial Biomarkers for Skin Cancer Immunotherapy

This project was made as a 10th grader and submitted to the 2020 IJAS science fair.

The dataset was taken from Matson et al's paper on the relative abundances of 63 bacteria from responder/non-responder data to anti-PD1 immunotherapy, found here: https://science.sciencemag.org/content/359/6371/104/tab-figures-data

Purpose: Melanoma is a type of skin cancer that makes up a large portion of skin cancer deaths. The development of immunotherapy––a treatment used to stimulate immune system function in cancer patients to destroy cancer cells––has been successfully applied in treating melanoma. However, many patients are still unresponsive to immunotherapy, thus understanding who may respond to treatment will be crucial. The purpose of my project is to identify key biomarkers in predicting the responders of immunotherapy patients from a published dataset of gut bacteria using machine learning techniques.

Procedure: Using a dataset of bacteria associated with immunotherapy outcomes, I created a customized random forest method to predict responders. I customized random forest––a machine learning approach based on multiple decision trees to achieve high prediction accuracy evaluated with the standard leave-one-out method. Specifically, I first sorted the bacteria by numerical importance in prediction, classified as Gini scores. I then applied random forest on each iteration set (the highest ranked bacteria recursively added per iteration), using leave-one-out cross validation sets. Thus, I obtained a subset of bacteria with the highest predictive accuracy.
