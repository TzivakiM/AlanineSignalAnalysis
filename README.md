# Alanine Signal Analysis

Sample files used for analysis of spectra produced by the BRuker EPR software. Input: ".asc" spectra files: B-field/Intensity measurement

Order of execution:

SignalHeights.py: Normalizes the x axis based on the B-field value of the standard and the x-axis based on the height of the standard intensity for each file in the containing folder
Output: csv file with signal heights (normalized as well as raw data)
errorasandplots.py: Performs the nomalization by packing density if provided with a suited input file. Can plot the data if uncommenting the relevant parts.
