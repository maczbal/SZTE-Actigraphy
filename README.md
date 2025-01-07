# SZTE-MIT-Actigraphy
Python-based web application to preprocess and convert triaxial acceleration signals into activity data, and to spectrally analyse them. The software is the work of Máté Miklós Perényi-Rimai for his bachelor thesis titled "Aktigráfiás jelek feldolgozását és spektrális analízisét végző webalkalmazás fejlesztése Python környezetben" (2024) under the supervision of Gergely Vadai and Bálint Maczák from the Department of Technical Infromatics, University of Szeged, Hungary.

# How to use

## Running the program
- Double-clicking on the executable file (" SZTE-MIT-Actigraphy .exe") will pop up a terminal window.
- The user interface will then automatically open as a web page in the default web browser. Closing the browser window does not stop the application, the user interface can be reopened from the browser via the web address "127.0.0.1:8055". The application can be shut down permanently by closing the terminal window.
## Loading a file
- On the web page, click on the "Select File" button to load the file to be scanned. The path (folder names) must not contain special characters.
    - Supported common file formats:
        - GENEActiv binary files based on the GGIRread (https://github.com/wadpac/GGIRread) R package, whose functionality was wrapped into a Python package (shared as installable in this repository).
        - EDF files through pyEDFlib Python package.
- Once the file is selected, the reading process starts immediately, indicated by a loading animation. The process typically takes 20-60 seconds. When the reading is finished, the raw acceleration data recorded along three axes is displayed. In the legend of the graph, you can see how the data of the different axes are represented in different colours, by clicking on one of its particular element you can hide/show the data of the given axis. To zoom in, select the magnifying glass icon in the toolbar in the top right corner, click in the graph and hold down the button while selecting the area to be zoomed in: dragging the cursor diagonally will zoom in on a rectangle, dragging horizontally or vertically will zoom in along the given axis. Click on the house icon to reset the view to default. For larger datasets, graphs may display a dynamically resampled data depending on the zoom level for optimal run speed, but this only affects visualisation, not any subsequent operations. By hovering the mouse pointer over a single data point, the exact time and amplitude value can be read.
## Preprocessing the acceleration data and optional activity calculation
- From the "Preprocessing Method" drop-down menu, you can choose which preprocessing you want to perform on the raw data before exporting it or calculating a Power Spectral Density (PSD) from it.
- Optional Activity Metric: From the "Activity Metric" drop-down menu, select an activity metric compatible with the preprocessing method specified above (e.g. ZCM metric for UFM preprocessing), which is used to calculate an activity value for a given length of time slices (epochs) of the preprocessed acceleration signal, thus obtaining an activity signal. The length of these epoch (the temporal resolution of the activity signal) can be specified in the input field "Epoch Length (s):". If you want to work with an acceleration signal, the activity metric must be left empty! For a detailed description of the different acceleration signal preprocessing techniques and activity metrics, see: https://doi.org/10.1371/journal.pone.0261718.
- By clicking on the "Plot" button, the signal set in the previous two steps is calculated and plotted.
## Exporting time series
- Use the "Start Date and Time (hour/min/sec)" and "End Date and Time (hour/min/sec)" time selectors to select the section you want to export or calculate a PSD from in the next step.
- Clicking on the "Export Time Series to CSV" button will export the time series and the associated metadata to a 1-1 file. In the file manager, an "Export" folder is created in the folder of the executable (SZTE-MIT-Actigraphy.exe), and the structured exported files are automatically included in this folder according to the file naming conventions.
## Spectral analysis and export of results
- By clicking on the "Spectrum" button, the PSD is calculated from the "Start Date and Time (hour/min/sec)" to "End Date and Time (hour/min/sec)" of the acceleration signal (or activity signal) generated in this way. At the same time, the fitting of the 1/f^β function (linear regression on the log transform of the filtered spectrum) is performed over the frequency range bounded by "Min Frequency Limit [Hz]" and "Max Frequency Limit [Hz]" of the averaged spectrum.
- If you want to change the fitting range, you can update the fitting by pressing the "Spectrum" button again after changing the frequency limits.
- The spectrum is shown on a logarithmic scale graph, where the original (not averaged) spectrum and the fitted line can be clicked in the legend. Above the graph you can also read the Slope (i.e., the exponent) and Y-Intercept of the fitted line and the accuracy of the fit as R^2.
- Similar to the export of time series, clicking on the "Export Spectrum to CSV" button will write the filtered spectrum and the fitted line, as well as the corresponding metadata, each to a file.
- Bins/decade defines the resolution of the log-binning-based averaging of the PSD (number of bins per decade), default is 20. For the detailed description of this averaging procedure, see Supplementary Material of https://doi.org/10.1038/s41598-024-52905-8.

# Description of exported metafiles and datafiles
## Naming scheme of the generated files
- For time series (TS)
    - TS_"unique_identifier"_ signaltype _cutstart_cutend_DT_META.csv
    - TS_"unique_identifier"_ signaltype _ cutstart _ cutend _dt_META.csv
- For PSDs (PSD)
    - PSD__"unique_identifier"_signaltype_ cutstart _ cutend _DT_FLO_FHI_BPD_META.csv
    - PSD__"unique_identifier"_ signaltype _ cutstart _ cutend _DT_FLO_FHI_BPD_DATA.csv
- Meaning of fields in filenames
    -  TS/PSD: Indicates whether it is a time series or a spectrum.
    -  META/DATA: Indicates whether it is a data file or a metafile.
    -  Unique identifier: a part extracted from the name or data of the input file. E.g., subject ID if given in the raw measurement data.
    -  Cutstart and Cutend: Indicates the time range that is selected on the acceleration or activity signal (start and end times, without year-month to shorten the file name). For TS, the data file contains the signal segment in this range, and for PSD, the spectrum (and fit) calculated from the signal segment in this range.
    -  DT: The time elapsed between successive data points in a temporal signal, this is 1/fs for an acceleration signal and the length of the epoch itself in seconds for an activity signal.
    -  FLO and FHI: for PSD, these are the lower and upper frequencies that define the range of the fitting, there is no such field for TS.
    -  BPD: Indicates the resolution of the bin-averaged spectrum (how many bins per decade), there is no such field for TS.
## Metafile structure (META)
### Fields both for TS and PSD
- File Name: the name of the processed raw acceleration data file.
- Measurement Sampling Rate [Hz]: the sampling frequency of the measurement.
- Acceleration Preprocessing Method: the method of preprocessing the acceleration signal (e.g. UFM).
- Activity Metric: the activity metric (e.g. ZCM) used to convert acceleration data to activity signal. For acceleration signal analysis this field is not relevant and empty.
- Activity Epoch Length [s]: the epoch length set for activity determination. For acceleration signal analysis this field is not relevant and empty.
- Examination Interval Start Time: examination interval start time (which is included in the filename but also with year and month here).
- Examination Interval End Time: similar to the previous one.
- Examination Interval Length: the total length of the cut-off period (hours, minutes, seconds).
### Extra fields for TS
- Timestamp of First Datapoint: if no data is recorded at the exact time of the start of the cut-off, the time stamp of the nearest subsequent datapoint is used.
- Elapsed Time Between Datapoints [s]: the time elapsed between successive data points in a time signal. Similar to the filename, for an acceleration signal 1/fs, and for an activity signal the length of the epoch itself in seconds.
### Extra fields for PSD
- Bins per Decade: the resolution of spectrum averaging.
- Fitting Interval LO Frequency [Hz]: the lower frequency of the fitting range.
- Fitting Interval HI Frequency [Hz]: similarly as above, the upper frequency of the fitting range.
- Slope: the slope of the fitted line (i.e. the estimated value of the exponent).
- Intercept: y-axis intercept of the fitted line.
- R2: Goodness of fit (R^2).
## Datafile structure (DATA)
- For TS, this csv file consists of 1 or 3 columns and contains only acceleration or activity values. For example, UFXYZ for preprocessing consists of 3 columns (1-1 column for each axis, order: x, y, z), but for UFM or activity signal only 1 column. A time column is not needed because it can be constructed from the timestamp of the first data point and the DT. These two information can be found in the associated metafile (META), see above.
- For PSD, the averaged spectrum is stored in 2 columns (one column is the frequency axis, the other is the magnitude values), and in addition, the frequency and magnitude values of the fitted line are stored in 2 additional columns.
