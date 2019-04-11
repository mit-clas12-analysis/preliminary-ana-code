# preliminary-ana-code
Code given to us by  Francois-Xavier Girod at the beginnning of April, 2019 to start on analysis

As per FX's email:

As we discussed during the meeting, please find attached some java code I wrote to analyze the channels of interests
As a disclaimer, I do not guarantee it to be without bug and I do not pretend that this is a model of well written code (it's not), however I hope this helps you start up

The argument passed are
- (1) a run number which is used to save the files in a directory called "plots[runNumber]
- (2) a text file containing on each line all the files you want to process, including directory path
- (3) an energy scale which sets the beam energy to 10.6 GeV for the scale > 8
- (4) the maximum number of events to process

The DVCS code actually has a couple more arguments, setting the number of bins in phi, and a flag for which calorimeter to use (0 for all, 1 for FTCal only, 2 for PCal only)

Please let me know if you have questions
I intend to start writing C++ / ROOT code next
