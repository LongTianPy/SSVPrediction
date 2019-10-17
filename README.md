# Splice Site Variation Prediction
Predicting whether a variation in Splice Site is pathogenic with an ensemble of alternative splicing site predictors.  
 
## Functions
- Given `chromosome` and `pos`, predict the scores of its probability of being a 
9-bp 5' splicing site and a 16-bp 3' splicing site by sliding 9-bp and 16-bp windows 
centered by this `pos`.
- Given `alt`, predict the scores of all potential splicing sites from the previous 
step. 
- Calculate the delta of scores.