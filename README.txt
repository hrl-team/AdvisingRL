README for the dataframe with raw data for Anllo & Hertz 2024

'data.frame':   34980 obs. of  32 variables:
 $ X.1           : int  1 2 3 4 5 6 7 8 9 10 ...                                     --> Placeholder
 $ id            : int  1 1 1 1 1 1 1 1 1 1 ...                                      --> Participant ID
 $ X             : int  5881 5882 5883 5884 5885 5886 5887 5888 5889 5890 ...        --> Participant ID (int. structure)   
 $ exptype       : chr  "No cost" "No cost" "No cost" "No cost" ...                  --> Advice consequence condition
 $ trialnum      : int  0 1 2 3 4 5 6 7 8 9 ...                                      --> Number of trial for the block
 $ choice        : int  0 1 1 1 1 1 1 1 0 0 ...                                      --> Correct (aka value-maximizing) 1, otherwise 0
 $ side          : int  -1 1 1 1 1 1 1 1 -1 -1 ...                                   --> -1 left, 1 right
 $ yesnochoice   : int  1 1 1 1 1 1 1 1 1 1 ...                                      --> Decision to disclose advice before seeing feedback (1 yes, 0 no)
 $ rt_yesnochoice: int  8811 1705 1313 1360 1113 1107 1163 1181 935 1755 ...         --> Advising decision response time in ms
 $ rt_choice     : int  5569 2676 1435 1266 1126 1609 1210 734 999 1402 ...          --> Choice response time in ms
 $ confidence    : int  -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 ...                            --> Deprecated column
 $ rt_confidence : int  -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 ...                            --> Deprecated column
 $ reward        : int  4 8 7 5 8 8 5 7 5 3 ...                                      --> Fish obtained
 $ adviserrate   : num  0 0 0 0 0 0 0 0 0 0 ...                                      --> Rating provided to adviser (condition "Reputation")
 $ blocknum      : int  0 0 0 0 0 0 0 0 0 0 ...                                      --> Block (by number)
 $ blocktype     : int  2 2 2 2 2 2 2 2 2 2 ...                                      --> Block (by difficulty: Easy, Medium, Hard)
 $ time          : int  194555 206910 215077 220885 225376 230335 235002 238951 243670 249087 ...  --> systime
 $ prev_reward   : int  0 4 8 7 5 8 8 5 7 5 ...                                      --> Previous reward 
 $ yesnochoicenum: int  1 1 1 1 1 1 1 1 1 1 ...                                      --> Decision to disclose advice before seeing feedback (1 yes, 0 no)
 $ PTM_all       : int  53 53 53 53 53 53 53 53 53 53 ...                            --> Score for PTM (all dimensions)
 $ PTM_anon      : int  11 11 11 11 11 11 11 11 11 11 ...                            --> Score for PTM (anon dimension)
 $ PTM_alt       : int  17 17 17 17 17 17 17 17 17 17 ...                            --> Score for PTM (altruism dimension)
 $ PTM_pub       : int  6 6 6 6 6 6 6 6 6 6 ...                                      --> Score for PTM (public dimension) 
 $ SPIN_all      : int  40 40 40 40 40 40 40 40 40 40 ...                            --> Score for SPIN (all dimensions)
 $ SPIN_fear     : int  19 19 19 19 19 19 19 19 19 19 ...                            --> Score for SPIN (fear dimension)
 $ SPIN_avoid    : int  14 14 14 14 14 14 14 14 14 14 ...                            --> Score for SPIN (avoidance dimension)
 $ SPIN_phys     : int  9 9 9 9 9 9 9 9 9 9 ...                                      --> Score for SPIN (phys dimension)
 $ age           : int  47 47 47 47 47 47 47 47 47 47 ...                            --> Participant age
 $ gender        : chr  "Male" "Male" "Male" "Male" ...                              --> Participant self-reported biological gender/sex
 $ scaledPTM     : num  -0.081 -0.081 -0.081 -0.081 -0.081 ...                       --> scaled PTM score (altruism dimension)
 $ scaledSPIN    : num  0.0333 0.0333 0.0333 0.0333 0.0333 ...                       --> scaled SPIN score (fear dimensions)
 $ id_2          : int  1 1 1 1 1 1 1 1 1 1 ...                                      --> Auxiliary ID
