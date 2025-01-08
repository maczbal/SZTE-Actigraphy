# SZTE-Actigraphy
Python-based web application to preprocess and convert actigraphic triaxial acceleration signals into activity data, and to spectrally analyse them. The software is the work of **Máté Miklós Perényi-Rimai** for his bachelor thesis titled "_Aktigráfiás jelek feldolgozását és spektrális analízisét végző webalkalmazás fejlesztése Python környezetben_" (2024) under the supervision of **Gergely Vadai** and **Bálint Maczák** from the **Department of Technical Infromatics, University of Szeged, Hungary**.

# Citing
If one uses this software, please refer to our article where we collected, categorized and compared the different activity calculation methods prevalent in the literature:<br/>**B. Maczák, G. Vadai, A. Dér, I. Szendi and Z. Gingl, "_Detailed analysis and comparison of different activity metrics_", PLOS ONE 16 (2021), e0261718, https://doi.org/10.1371/journal.pone.0261718**

# How to use
The program can be converted into stand-alone executable with PyInstaller by executing the following command in terminal:<br/>
```python -m PyInstaller --onefile --name “SZTE-MIT-Actigraphy” "Actigraphy.py"```
## Running the program
- Double-clicking on the created executable file (i.e., _SZTE-MIT-Actigraphy.exe_) will pop up a terminal window.
- The user interface will then automatically open as a web page in the default web browser. Closing the browser window does not stop the application, the user interface can be reopened from the browser via the web address "127.0.0.1:8055". The application can be shut down permanently by closing the terminal window.
## Loading a file
- On the web page, click on the _**Select File**_ button to select the file to be processed. The path of the file must not contain special characters.
    - Supported common file formats:
        - GENEActiv binary files based on the GGIRread (https://github.com/wadpac/GGIRread) R package, whose functionality was wrapped into a Python package and shared in the folder _GENEActiv Reader Package_. **You must install this package before creating the executable**, for this, navigate into the folder and execute the following command in terminal: <br/>```pip install .```
        - EDF files through pyEDFlib Python package.
- Once the file is selected, the reading process starts immediately, indicated by a loading animation. The process typically takes 20-60 seconds. When the reading is finished, the raw acceleration data recorded along three axes is displayed. In the legend of the graph, you can see how the data of the different axes are represented in different colours, by clicking on one of its particular element you can hide/show the data of the given axis. To zoom in, select the magnifying glass icon in the toolbar of the graph in the top right corner, click in the graph and hold down the button while selecting the area to be zoomed in: dragging the cursor diagonally will zoom in on a rectangle, dragging horizontally or vertically will zoom in along the given axis. Click on the house icon to reset the view to default. For larger datasets, graphs may display a dynamically resampled data depending on the zoom level for optimal run speed, but this only affects visualisation, not any subsequent operations. By hovering the mouse pointer over a single data point, the exact time and amplitude value can be read.
## Preprocessing the acceleration data and optional activity calculation
- From the _**Preprocessing Method**_ drop-down menu, you can choose which preprocessing you want to perform on the raw data before exporting it or calculating a Power Spectral Density (PSD) from it.
- Optionally, from the _**Activity Metric**_ drop-down menu, select an activity metric compatible with the preprocessing method specified above (e.g. ZCM metric for UFM preprocessing), which is used to calculate an activity value for a given length of consecutive, non-overlapping time slices (epochs) of the preprocessed acceleration signal, thus obtaining an activity signal. The length of these epoch (the temporal resolution of the activity signal) can be specified in the input field _**Epoch Length (s)**_. If you want to work with an acceleration signal, the activity metric must be left empty! For a detailed description of the different acceleration signal preprocessing techniques and activity metrics, and their compatibility, see our previous work [[1]](#1).
- By clicking on the _**Plot**_ button, the signal defined in the previous two steps is calculated and plotted.
## Exporting time series
- Use the _**Start Date and Time (hour/min/sec)**_ and _**End Date and Time (hour/min/sec)**_ time selectors to select the interval you want to export or calculate a PSD from in the next step.
- Clicking on the _**Export Time Series to CSV**_ button will export the time series and the associated metadata to a 1-1 file (the description of these files are given below). Consequently, an _**Export**_ folder is created next to the executable, and the structured exported files are automatically generated in this folder according to the file naming conventions (also described below).
## Spectral analysis and export of results
- By clicking on the _**Spectrum**_ button, the PSD is calculated from the data spanning from _**Start Date and Time (hour/min/sec)**_ to _**End Date and Time (hour/min/sec)**_ of the acceleration (or activity) signal. At the same time, the fitting of the $S(f) \propto f^{-\beta}$ function (using linear regression over the log transform spectrum) is performed over the frequency range bounded by _**Min Frequency Limit [Hz]**_ and _**Max Frequency Limit [Hz]**_ of the _averaged_ PSD. For the detailed description of this log-binning-based spectrum averaging procedure, see the Supplementary Material of our previous work [[2]](#2).
- If you want to change the fitting range, you can update the fitting by pressing the _**Spectrum**_ button again after changing the frequency limits.
- The spectrum is shown on a logarithmic scale graph, where the original (not averaged) spectrum and the fitted line can be clicked in the legend to make them visible.
- Above the graph you can also read the _**Slope**_ (i.e., the $\beta$ exponent) and _**Y-Intercept**_ of the fitted line and the accuracy of the fit as $R^2$.
- Similar to the export of time series, clicking on the _**Export Spectrum to CSV**_ button will write the averaged spectrum and the fitted line, as well as the corresponding metadata, each to a file.
- _**Bins/decade**_ defines the resolution of the log-binning-based averaging of the PSD (number of bins per decade), default is 20.

# Description of exported metafiles and datafiles
## Naming scheme of the generated files
- For time series (TS)
    - TS_unique_identifier_signaltype_cutstart_cutend_DT_META.csv
    - TS_unique_identifier_signaltype_cutstart_cutend_dt_META.csv
- For PSDs (PSD)
    - PSD_unique_identifier_signaltype_cutstart_cutend_DT_FLO_FHI_BPD_META.csv
    - PSD_unique_identifier_ signaltype_cutstart_cutend_DT_FLO_FHI_BPD_DATA.csv
- Meaning of fields in filenames
    -  _**TS**_/_**PSD**_: Indicates whether it is a time series or a spectrum.
    -  _**META**_/_**DATA**_: Indicates whether it is a data file or a metafile.
    -  _**Unique identifier**_: A part extracted from the name or data of the input file. E.g., subject ID if given in the raw measurement data.
    -  _**Cutstart**_ and _**Cutend**_: Indicates the time range that is selected on the acceleration or activity signal (start and end times, without year-month to shorten the file name). For TS, the data file contains the signal segment in this range, and for PSD, the spectrum (and fit) calculated from the signal segment in this range.
    -  _**DT**_: The time elapsed between successive data points in a temporal signal, this is $1\fs$ for an acceleration signal and the length of the epoch itself in seconds for an activity signal.
    -  _**FLO**_ and _**FHI**_: For PSD, these are the lower and upper frequencies that define the range of the fitting, there is no such field for TS.
    -  _**BPD**_: Indicates the resolution of the bin-averaged spectrum (how many bins per decade), there is no such field for TS.
## Metafile structure (META)
### Fields both for TS and PSD
- _**File Name**_: The name of the processed raw acceleration data file.
- _**Measurement Sampling Rate [Hz]**_: The sampling frequency of the measurement.
- _**Acceleration Preprocessing Method**_: The method of preprocessing the acceleration signal (e.g. UFM).
- _**Activity Metric**_: The activity metric (e.g. ZCM) used to convert acceleration data to activity signal. For acceleration signal analysis this field is not relevant and empty.
- _**Activity Epoch Length [s]**_: The epoch length set for activity determination. For acceleration signal analysis this field is not relevant and empty.
- _**Examination Interval Start Time**_: Examination interval start time (which is included in the filename but also with year and month here).
- _**Examination Interval End Time**_: Similar to the previous one.
- _**Examination Interval Length**_: The total length of the cut-off period (hours, minutes, seconds).
### Extra fields for TS
- _**Timestamp of First Datapoint**_: If no data is recorded at the exact time of the start of the cut-off, the time stamp of the nearest subsequent datapoint is used.
- _**Elapsed Time Between Datapoints [s]**_: The time elapsed between successive data points in a time signal. Similar to the filename, for an acceleration signal $1\fs$, and for an activity signal the length of the epoch itself in seconds.
### Extra fields for PSD
- _**Bins per Decade**_: The resolution of spectrum averaging.
- _**Fitting Interval LO Frequency [Hz]**_: The lower frequency of the fitting range.
- _**Fitting Interval HI Frequency [Hz]**_: Similarly as above, the upper frequency of the fitting range.
- _**Slope**_: The slope of the fitted line (i.e. the estimated value of the exponent).
- _**Intercept**_: Y-axis intercept of the fitted line.
- _**R2**_: Goodness of fit ($R^2$).
## Datafile structure (DATA)
- For TS, this csv file consists of 1 or 3 columns and contains only acceleration or activity values. For example, exported acceleration data preprocessed using UFXYZ consists of 3 columns (1-1 column for each axis, order: x, y, z), but for UFM or activity signals, it is only 1 column of amplitude values. A time column is not needed because it can be constructed from the timestamp of the first data point and the DT. These two information can be found in the associated metafile (META), see above.
- For PSD, the averaged spectrum is stored in 2 columns (one column is the frequency axis, the other is the magnitude values), and in addition, the frequency and magnitude values of the fitted line are stored in 2 additional columns.

## References
<a id="1">[1]</a> 
B. Maczák, G. Vadai, A. Dér, I. Szendi and Z. Gingl, "_Detailed analysis and comparison of different activity metrics_", PLOS ONE 16 (2021), e0261718, https://doi.org/10.1371/journal.pone.0261718 <br/>
<a id="2">[2]</a> 
B. Maczák, Z. Gingl and G. Vadai, "_General spectral characteristics of human activity and its inherent scale-free fluctuations_", SCIENTIFIC REPORTS 14 (2024), 2604, https://doi.org/10.1038/s41598-024-52905-8
