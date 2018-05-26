# Bootstrap Magnitude Based Predictions (BMBP)
Copyright 2018 Mladen Jovanovic

For videos visit: http://www.complementarytraining.com/bmbp

There has been a lot of fuss regarding Magnitude Based Inference (MBI) lately, peaking with the recent article and video by Kristin Sainani. The originators of MBI, Will Hopkins and Alan Batterham, recently posted their commentary of the critiques at sportsci.org. With the video below I wanted to give my opinion, but also provide rationale for “novel” techniques, that might seem more intuitive for the coaches. I am aware that I might end up labelled stupid or erroneous, but that is the risk I am willing to take.

Long story short, Null Hypothesis Testing (NHT) and Magnitude Based Inference (MBI) revolve around estimating MEAN/AVERAGE effects in population. But, as a coach I am not really concerned about estimates of the mean/average effects in the populations, but INDIVIDUAL effects in my group (and/or population). So, assuming we know the Smallest Worthwhile Change (SWC), which could be estimated using 0.2 x SD between individuals, or using our best guess as a coach (e.g. improvement in 1RM of 2.5kg might be deem smallest worthwhile improvement, while 1kg is not), we might be interested in estimating how many ATHLETES are likely to get harmful/trivial/beneficial effects. Not mean/average effects, but individual effects.

So, to solve this issues, I leaned more toward “predictive camp” of statistical modeling and used very useful (but computationally intensive) technique called “bootstrap”. In the video below I am explaining (and providing free R code) how to get a glimpse if a certain test/measurement is reliable enough by using SWC-to-Typical Error (TE) metric, as well as a “novel” metric that I’ve called “HitRate”, which are estimates of signal vs noise. After estimating TE, we can proceed in checking effects of the intervention. Here I provide a bootstrap estimates of both MBIs, as well as “novel” bootstrap magnitude based predictions (BMBP), and also provide a way how to utilize TE in the bootstrap that can give us more realistic estimates of the effects.

I believe that Magnitude Based Predictions answer practical questions from the field (which are more predictive rather than inferential) much better than Magnitude Based Inferences. For example, coach might asks “What is the probability that Jorge will improve with this intervention?”. With BMBP you have a tool to answer that question. Because in real life, practitioners are interested in N=1, not mean effects in the population.

If you plan using these techniques or referencing this article, please use the following reference:

Jovanovic, M.(2018). Bootstrap Magnitude Based Predictions. Complementary Training. Retrieved from http://www.complementarytraining.com/bmbp
