---
title: "Quiz 2"
author: "Monikrishna Roy"
date: "9/18/2021"
output: 
  pdf_document:
    keep_md: yes
---



## Question 1
Download GSE9782.RData from UNR Canvas. Use the paste0() function to provide the complete path of the file stored in your working folder. Use the load() function to load the data into your RStudio working environment.


```r
load(paste0("Project-1/GSE9782.RData"))
```

## Question 2
Use the function which() to determine the indices of patients who were treated, and the indices of patients who were not treated. Note that the treatment information is stored in the group object.


```r
treated <- which(group == "treated")
untreated <- which(group != "treated")
```

## Question 3
Compute the following matrices:

a) log base 2 of the expression data
b) log base 10 of the expression data
c) log base 5 of the expression data


```r
log_2_data <- log(x = dataGSE9782, base = 2)
log_10_data <- log(x = dataGSE9782, base = 10)
log_5_data <- log(x = dataGSE9782, base = 5)
```

## Question 4
Compute the difference in mean log base 2 expression between the two conditions (treated versus control) using the for loop. Name the variable of the results as logFC1


```r
# a function to calculate the difference in mean log base 2
cal_mean_diff <- function(variables) {
mean(log2(variables[treated])) - mean(log2(variables[untreated]))
}
# used for to call the function
logFC1 = vector("numeric", nrow(dataGSE9782))
for(r in 1:nrow(dataGSE9782)) {
  logFC1[r] = cal_mean_diff(dataGSE9782[r, ])
}
names(logFC1) <- rownames(dataGSE9782)
```

## Question 5
Compute the difference in mean log base 2 expression between the two conditions (treated versus control) using the apply() function. Name the variable of the results as logFC2. Make sure that logFC1 and logFC2 have exactly the same values.


```r
# a function to calculate the difference in mean log base 2
cal_mean_diff <- function(variables) {
mean(log2(variables[treated])) - mean(log2(variables[untreated]))
}
# used apply to call the function
logFC2 <- apply(dataGSE9782, 1, cal_mean_diff)

sum(logFC2 != logFC1)
```

```
## [1] 0
```
