# %%
from datetime import datetime, timezone, timedelta
from plotly_resampler import FigureWidgetResampler
import plotly.graph_objects as go
from scipy import stats, signal
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd

import base64
import os
from dash import Dash, Input, Output, callback_context, dcc, html, no_update, State
from dash.dependencies import Input, Output
import pandas as pd
import numpy as np
import plotly.graph_objs as go
from plotly_resampler import FigureResampler
from plotly.subplots import make_subplots
from tkinter import Tk, filedialog
from pybindRetTest import GENEActivReader
from pyedflib import highlevel
import gc
import re

import webbrowser
from threading import Timer
from MIT_reader import actigraph_read_binary_copy


# %%
# Initialize Dash application and declare global variables
app = Dash(__name__, suppress_callback_exceptions=True)

# Resampler figures for data visualization
fig_xyz: FigureResampler = FigureResampler()
fig: FigureResampler = FigureResampler()
fig_spectrum : FigureResampler = FigureResampler()

# Global variables for data processing and configurations
bins_per_decade_global = 20
epoch_length = 60  
global last_df, actual_df, file_type, patientcode, UFNM_df, UFM_df, TAT_deltaT
last_df = actual_df = file_type = patientcode = UFNM_df = UFM_df = TAT_deltaT = None

global FMpre_df, FMpost_df, FXYZ_df, FX_df, FY_df, FZ_df
FMpre_df = FMpost_df = FXYZ_df = FX_df = FY_df = FZ_df = None

global FX_df_mag, FY_df_mag, FZ_df_mag, UFX_df, UFY_df, UFZ_df
FX_df_mag = FY_df_mag = FZ_df_mag = UFX_df = UFY_df = UFZ_df = None

global UFX_df_mag, UFY_df_mag, UFZ_df_mag, HFMpre_df, measurement_freq, deltaT
UFX_df_mag = UFY_df_mag = UFZ_df_mag = HFMpre_df = measurement_freq = deltaT = None




# %%
# Function to create a DataFrame from raw acceleration data
def df(x_acc, y_acc, z_acc, start_time, deltaT):
    measurement_freq = round(deltaT.total_seconds()*1000, 13)
    mag_acc = np.sqrt(x_acc**2 + y_acc**2 + z_acc**2)
    timestamps = pd.date_range(start_time, periods=x_acc.shape[0], freq=(str(measurement_freq) + "ms"))
    timestamps_without_tz = timestamps.tz_localize(None)
    df = pd.DataFrame({'Timestamp': timestamps_without_tz.values, 'X-Axis': x_acc, 'Y-Axis': y_acc, 'Z-Axis': z_acc, 'Magnitude': mag_acc})
    return df


# %%
# Function to create epochs from a DataFrame
def epoch(df, epoch_s):
    start_time = pd.Timestamp(df['Timestamp'].iloc[0]).ceil('min')
    end_time = pd.Timestamp(df['Timestamp'].iloc[-1])
    mintervals = pd.date_range(start=start_time, end=end_time, freq=str(epoch_s) + "s")
    epoch_df = pd.DataFrame({'Timestamp': df['Timestamp'], 'Magnitude': df['Magnitude']})
    epoch_df['Timestamp'] = pd.cut(epoch_df['Timestamp'], bins=mintervals, right=False, labels=mintervals[:-1]).astype('datetime64[ns]')
    epoch_df = epoch_df.dropna()
    return epoch_df


# %%
# Functions of the different activity metrics

def ZCM(data, threshold):
    signs = np.sign(data - threshold)  
    zero_crossings = np.count_nonzero(np.diff(signs) != 0)  
    return zero_crossings

# %%
def TAT(data, threshold, deltaTime):
    num_above_threshold = np.count_nonzero(data > threshold)*deltaTime
    return num_above_threshold

# %%
def MAD(data):
    mean_value = np.mean(data)
    mad_value = np.mean(np.abs(data - mean_value))
    return mad_value

# %%
def ENMO(data):
    enmo_value = np.mean(np.maximum(data - 1, 0))
    return enmo_value

# %%
def PIM(data, epoch_length, delta_time, method='-'):
   
    epoch_df = epoch(data, epoch_length)
    
    if method in ['UFM']:
        raw_pim_values = delta_time * epoch_df['Magnitude'].groupby(epoch_df['Timestamp']).sum()
        gravity_correction = delta_time * epoch_df['Magnitude'].groupby(epoch_df['Timestamp']).count()
        
        corrected_pim_values = np.abs(raw_pim_values - gravity_correction)
    elif method in ['FX','FZ','FY', 'FMpost']:
        corrected_magnitude = np.abs(epoch_df['Magnitude'])
        corrected_pim_values = delta_time * corrected_magnitude.groupby(epoch_df['Timestamp']).sum()
    else:
        corrected_magnitude = epoch_df['Magnitude']
        corrected_pim_values = delta_time * corrected_magnitude.groupby(epoch_df['Timestamp']).sum()

    pim_df = pd.DataFrame({
        'Timestamp': corrected_pim_values.index,
        'Activity': corrected_pim_values.values
    }).reset_index(drop=True)
    return pim_df

# %%
def HFEN(data):
    hfen_value = np.mean(data)
    return hfen_value

# %%
# Power Spectral Density (PSD) calculation
def psd(t,x):
    N = len(x)
    dt = (t[1] - t[0]).total_seconds()
    fs = 1/dt
    T = dt*N
    df = 1/T
    a = (np.abs(np.fft.fft(x.to_numpy()))**2) / (N**2)
    iq = int(np.floor(N/2))
    r = N-2*iq
    a = a[0:r+iq]
    a[1:r+iq] *= 2
    a *= T
    f = np.arange(start = 0, stop = fs/2, step=df)
    return pd.DataFrame({'Frequency':f, 'Magnitude':a})

# %%
# Fit a linear regression to a section of the PSD
def fit_line_to_psd(df_psd_reduced, fit_frequency_bound_lo=0.0001, fit_frequency_bound_hi=0.01):

    # Filter the PSD data for the specified frequency range
    df_psd_to_fit = df_psd_reduced[(df_psd_reduced['Frequency'] > fit_frequency_bound_lo) & (df_psd_reduced['Frequency'] <= fit_frequency_bound_hi)]

    # Perform a linear regression on log-transformed data
    linregress_res = stats.linregress(np.log10(df_psd_to_fit['Frequency']), np.log10(df_psd_to_fit['Magnitude']))

    # Calculate the fitted line
    df_fit = pd.DataFrame({'Frequency': df_psd_to_fit['Frequency'], 'Magnitude': np.power(10, linregress_res.slope * np.log10(df_psd_to_fit['Frequency']) + linregress_res.intercept)})
    return df_fit, linregress_res

# %%
def butter_bandpass_filter(data, fs):
    b, a = signal.iirfilter(3, Wn=[0.25, 2.5], btype="bandpass", ftype="butter", fs=fs)
    y = signal.lfilter(b, a, data)
    return y

# %%
def highpass_filter(data, cutoff, fs, order=4):
    nyquist = 0.5 * fs
    normal_cutoff = cutoff / nyquist
    b, a = signal.butter(order, normal_cutoff, btype='high', analog=False)
    filtered_data = signal.lfilter(b, a, data)
    return filtered_data

# %%
# Reset global variables to their initial state, used when a new file is selected
def reset_global_variables():
    global last_df, actual_df, UFM_df, UFNM_df, FMpre_df, FMpost_df, FXYZ_df
    global FX_df, FY_df, FZ_df, FX_df_mag, FY_df_mag, FZ_df_mag
    global UFX_df, UFY_df, UFZ_df, UFX_df_mag, UFY_df_mag, UFZ_df_mag
    global measurement_freq, deltaT, deltaTime, file_type, patientcode
    global fig_xyz, fig, fig_spectrum

    last_df = actual_df = UFM_df = UFNM_df = FMpre_df = FMpost_df = FXYZ_df = FX_df = FY_df = FZ_df = FX_df_mag = FY_df_mag = FZ_df_mag = UFX_df = UFY_df = UFZ_df = UFX_df_mag = UFY_df_mag = UFZ_df_mag = measurement_freq = deltaT = deltaTime = file_type = patientcode = None

    fig_xyz.replace(go.Figure())
    fig.replace(go.Figure())
    fig_spectrum.replace(go.Figure())
    gc.collect()  


# %%
# Dropdown options for preprocessing methods
first_dropdown_options = [
    {'label': 'UFXYZ', 'value': 'UFXYZ'},
    {'label': 'FXYZ', 'value': 'FXYZ'},
    {'label': 'UFM', 'value': 'UFM'},
    {'label': 'UFNM', 'value': 'UFNM'},
    {'label': 'FMpost', 'value': 'FMpost'},
    {'label': 'FMpre', 'value': 'FMpre'},
    {'label': 'HFMpre', 'value': 'HFMpre'},
    {'label': 'UFX', 'value': 'UFX'},
    {'label': 'UFY', 'value': 'UFY'},
    {'label': 'UFZ', 'value': 'UFZ'},
    {'label': 'FX', 'value': 'FX'},
    {'label': 'FY', 'value': 'FY'},
    {'label': 'FZ', 'value': 'FZ'},
    
]

# Dropdown options for activity metrics
second_dropdown_options = [
    {'label': 'ZCM', 'value': 'ZCM'},
    {'label': 'TAT', 'value': 'TAT'},
    {'label': 'MAD', 'value': 'MAD'},
    {'label': 'ENMO', 'value': 'ENMO'},  
    {'label': 'PIM', 'value': 'PIM'},
    {'label': 'HFEN', 'value': 'HFEN'},
    {'label': '-', 'value': '-'}
]

# Layout configuration for the application
app.layout = html.Div([
    html.Div([
        html.Div([
            html.Label('File Path: ', style={'color': '#ffffff', 'font-size': '22px'}),
        ]),
        
        html.Div([
            
            dcc.Input(id='file-path-input', type='text', placeholder='Select a file', style={'width': '20%', 'height': '30px', 'font-size': '18px', 'marginBottom': '15px'}),
            html.Button('Select File', id='select-file-button', n_clicks=0, style={'marginLeft': '10px', 'backgroundColor': '#4CAF50', 'color': '#ffffff', 'padding': '10px 15px', 'font-size': '20px', 'marginBottom': '15px'}),
            html.Label('EDF Channels:', style={'color': '#ffffff', 'font-size': '25px', 'marginLeft': '30px' }),
            html.Div([
            
            dcc.Checklist(
                id='channel-checklist',
                options=[
                    {'label': 'CH1', 'value': 0},
                    {'label': 'CH2', 'value': 1},
                    {'label': 'CH3', 'value': 2},
                    {'label': 'CH4', 'value': 3},
                    {'label': 'CH5', 'value': 4},
                    {'label': 'CH6', 'value': 5}
                ],
                value=[3,4,5],  
                inline=True,
                labelStyle={'display': 'inline-block', 'marginRight': '10px', 'fontSize': '18px'},
                style={ 'color': '#ffffff', 'marginLeft' : '10px', 'marginTop':'12px' ,'fontSize': '18px'}, 
            )
        ], style={'display': 'inline-block', 'verticalAlign': 'middle', 'marginBottom': '15px'}),
            dcc.Loading(id="loading-spinner-xyz1", type="circle", className="custom-loading-spinner", children=[
            html.Div(id="graph-container-xyz")
        ], style={'display': 'flex', 'justifyContent': 'left', 'alignItems': 'center', 'margin': '10px 0'}),
        
        html.Div(id='output-data-upload', style={'textAlign': 'center', 'margin': '10px 0', 'color': '#ffffff'}),

        ]),

        

        html.Hr(style={'borderColor': '#ffffff'}),
        
        
        
        html.Div([
            html.Div([
                html.Label('Preprocessing Method:', style={'color': '#ffffff', 'font-size': '20px', }),
                dcc.Dropdown(
                    id='preprocessing-method',
                    options=first_dropdown_options,
                    value=None,
                    searchable=False,
                    clearable=False,
                    style={'width': '100%', 'marginTop' : '5px', 'marginBottom' : '5px'}
                ),
                html.Label('Activity Metric:', style={'color': '#ffffff', 'font-size': '20px',}),
                dcc.Dropdown(
                    id='activity-method',
                    options=second_dropdown_options,
                    value='-',
                    disabled=True,
                    searchable=False,
                    clearable=False,
                    style={'width': '100%', 'marginTop' : '5px'}
                ),
                html.Label('Epoch Length (s):', style={'color': '#ffffff', 'font-size': '20px'}),
                dcc.Input(id='epoch-length-input', type='number', placeholder='Sec', value=60, style={'width': '20%', 'height': '22px', 'font-size': '20px'}),
                html.Button("Plot", id="plot-button", n_clicks=0, style={'marginRight': '10px', 'backgroundColor': '#4CAF50', 'color': '#ffffff', 'padding': '10px 30px', 'font-size': '20px' ,'marginTop' : '10px', 'marginLeft' : '10px', 'marginBottom': '20px'}),
                



        html.Hr(style={'borderColor': '#ffffff', 'marginBottom': '20px'}),


        html.Label('Export Time Series', style={'color': '#ffffff', 'font-size': '30px', 'marginBottom' : '30px', 'padding': '42px',  }),


        html.Div(style={'margin-top': '20px'}, children=[
            html.Label('Start Date and Time (hour/min/sec):', style={'color': '#ffffff', 'font-size': '20px', 'marginTop' : '40px' }),
        ]),
        html.Div([
            dcc.DatePickerSingle(
                id='start-date-picker',
                display_format='YYYY-MM-DD',
                persisted_props=['date'],
                persistence_type='session',
                style={'marginRight': '5px', 'marginTop' : '5px' }
                ),
            dcc.Input(id='start-hour', type='number', placeholder='Hour', min=0, max=23, value=0, style={'width': '50px', 'marginRight': '5px' , 'height' : '25px'}),
            dcc.Input(id='start-minute', type='number', placeholder='Minute', min=0, max=59, value=0, style={'width': '50px', 'marginRight': '5px', 'height' : '25px'}),
            dcc.Input(id='start-second', type='number', placeholder='Second', min=0, max=59, value=0, style={'width': '50px' , 'height' : '25px'}),
        ]),
        html.Div([
            html.Label('End Date and Time (hour/min/sec):', style={'color': '#ffffff', 'font-size': '20px'}),
        ]),
        html.Div([
            dcc.DatePickerSingle(
                id='end-date-picker',
                display_format='YYYY-MM-DD',
                persisted_props=['date'],
                persistence_type='session',
                style={'marginRight': '5px', 'marginTop' : '5px' }
                ),
            dcc.Input(id='end-hour', type='number', placeholder='Hour', min=0, max=23, value=0, style={'width': '50px', 'marginRight': '5px' , 'height' : '25px'}),
            dcc.Input(id='end-minute', type='number', placeholder='Minute', min=0, max=59, value=0, style={'width': '50px', 'marginRight': '5px' , 'height' : '25px'}),
            dcc.Input(id='end-second', type='number', placeholder='Second', min=0, max=59, value=0, style={'width': '50px' , 'height' : '25px'}),

            html.Div([
            html.Button("Export Time Series to CSV", id="export-csv-button", n_clicks=0, style={'backgroundColor': '#4CAF50', 'color': '#ffffff', 'padding': '10px 20px', 'font-size': '20px', 'width' : '300px'}), 
            ], style={  'padding': '40px'}),
        ]),
                
            ], style={'width': '20%', 'marginRight': '15px', 'padding': '20px'}),  
            html.Div([
                dcc.Loading(id="loading-spinner-xyz2", type="circle", className="custom-loading-spinner", children=[
                    dcc.Graph(id="graph-id", style={'height': '600px'})  
                ]),
            ], style={'width': '80%', 'padding': '20px'})
        ], style={'display': 'flex'}),

       

        html.A(id="download-link", download="activity_data.csv", href="", target="_blank", style={'display': 'block', 'textAlign': 'center', 'margin': '10px 0', 'color': '#ffffff'}),
        html.A(id="download-link-reduced", download="psd_reduced.csv", href="", target="_blank", style={'display': 'block', 'textAlign': 'center', 'margin': '10px 0', 'color': '#ffffff'}),
        html.Hr(style={'borderColor': '#ffffff'}),
        
        html.Div([
            html.Button("Spectrum", id="calculate-spectrum-button", n_clicks=0, style={'backgroundColor': '#4CAF50', 'color': '#ffffff', 'padding': '10px 20px', 'font-size': '20px', 'marginRight':'15px'}),
            html.Button("Export Spectrum to CSV", id="export-csv-reduced-button", n_clicks=0, style={'backgroundColor': '#4CAF50', 'color': '#ffffff', 'padding': '10px 20px', 'font-size': '20px', 'marginRight' : '10px'}), 
            html.Label('Fit:', style={'color': '#ffffff', 'font-size': '26px', 'marginTop': '10px', 'marginRight': '15px'}),
            html.Label('Min Frequency Limit [Hz]:', style={'color': '#ffffff', 'font-size': '20px', 'marginTop': '10px'}),
            dcc.Input(id='fit-frequency-bound-lo', type='number', placeholder='Fit Start', value=0.0001, style={'width': '60px', 'marginRight': '10px' , 'height' : '25px'}),
            html.Label('Max Frequency Limit [Hz]:', style={'color': '#ffffff', 'font-size': '20px', 'marginTop': '10px'}),
            dcc.Input(id='fit-frequency-bound-hi', type='number', placeholder='Fit End', value=0.01, style={'width': '60px', 'marginRight': '5px' , 'height' : '25px'}),
            html.Label('Bins/decade:', style={'color': '#ffffff', 'font-size': '20px', 'marginTop': '10px'}),
            dcc.Input(id='bins-per-decade', type='number', placeholder='', value=20, style={'width': '60px', 'marginRight': '10px' , 'height' : '25px'}),

        ], style={'textAlign': 'left', 'margin': '10px 0'}),

        html.Div([
            
            html.Label('Slope:', style={'color': '#ffffff', 'font-size': '20px', 'marginTop': '10px'}),
            dcc.Input(id='slope', type='text', readOnly=True, style={'width': '60px', 'marginRight': '13px' , 'height' : '25px'}),
            html.Label('Y-Intercept:', style={'color': '#ffffff', 'font-size': '20px', 'marginTop': '10px'}),
            dcc.Input(id='intercept', type='text', readOnly=True, style={'width': '145px', 'marginRight': '13px' , 'height' : '25px'}),
            html.Label('R²:', style={'color': '#ffffff', 'font-size': '20px', 'marginTop': '10px'}),
            dcc.Input(id='r2', type='text',  readOnly=True, style={'width': '145px', 'marginRight': '5px' , 'height' : '25px'}),
            
        ], style={'textAlign': 'left', 'margin': '10px 0'}),
        
        html.Div(id='spectrum-graph-container', children=[
            dcc.Loading(id="loading-spinner-xyz3", type="circle", className="custom-loading-spinner", children=[
                dcc.Graph(id='spectrum-graph')
            ]),
        ], style={'margin': '10px 0'}),
        
    ], style={'width': '97%', 'margin': '0 auto', 'padding': '20px', 'backgroundColor': '#333333', 'borderRadius': '10px'})
])

# %%
# Directory for the exported data
DOWNLOAD_DIR = "Export"

# Create the directory if it doesn't exist
if not os.path.exists(DOWNLOAD_DIR):
    os.makedirs(DOWNLOAD_DIR)

# Function to open a file selection dialog using Tkinter
def select_file():
    file_path = None
    root = Tk()
    root.withdraw()  
    root.call('wm', 'attributes', '.', '-topmost', True)  
    file_path = filedialog.askopenfilename()  
    root.destroy()  
    return file_path


@app.callback(
    Output('file-path-input', 'value'),
    [Input('select-file-button', 'n_clicks')]
)
def update_file_path_input(n_clicks):
    if n_clicks > 0:
        global file_path
        file_path = select_file()  
        return file_path  


# Function to update the first graph when a file is uploaded
# %%
@app.callback(
    Output('output-data-upload', 'children'),
    Output('graph-container-xyz', 'children'), 
    [Input('file-path-input', 'value'),
     State('channel-checklist', 'value')]
)
def update_graph(file_path, selected_channels):

    if not file_path:
        return no_update, no_update
    
    reset_global_variables()
    global UFM_df,csv_file_name, measurement_freq, deltaT, deltaTime, file_type

    csv_file_name = file_path.split('/')
    string = file_path.split('.')

    if file_path and string[-1] == "actigraph":
        file_type = "actigraph"
        x_acc, y_acc, z_acc, start_time, deltaT = actigraph_read_binary_copy(file_path, 'Europe/Budapest')
        deltaTime = deltaT.total_seconds()

        UFM_df = df(x_acc, y_acc, z_acc, start_time, deltaT)

        del x_acc, y_acc, z_acc
        gc.collect()
        
        fig_xyz.replace(go.Figure())
        fig_xyz.add_trace(go.Scatter(showlegend=True, name='X'), hf_x=UFM_df['Timestamp'], hf_y=UFM_df['X-Axis'], max_n_samples=5000)
        fig_xyz.add_trace(go.Scatter(showlegend=True, name='Y'), hf_x=UFM_df['Timestamp'], hf_y=UFM_df['Y-Axis'], max_n_samples=5000)
        fig_xyz.add_trace(go.Scatter(showlegend=True, name='Z'), hf_x=UFM_df['Timestamp'], hf_y=UFM_df['Z-Axis'], max_n_samples=5000)
        fig_xyz.update_layout(height=600, xaxis_title="Time", yaxis_title='Acceleration [g]', margin=dict(l=37, r=20, t=10, b=10 ), legend=dict(orientation="h"))
        
        graph_component = dcc.Graph(figure=fig_xyz, id="xyz")

        return "", graph_component
    
    elif file_path and string[-1] == "bin":
        
        file_type = "bin"
        file_path = file_path.replace('\\','\\\\')
        ret = GENEActivReader(file_path)
        x_acc = ret.getAccelerationX()
        y_acc = ret.getAccelerationY()
        z_acc = ret.getAccelerationZ()
        mag_acc = np.sqrt(x_acc**2 + y_acc**2 + z_acc**2)
        start_time = ret.getMeasurementStartTime()
        parsed_start_time = pd.Timestamp(datetime.strptime(start_time, "%Y-%m-%d %H:%M:%S.%f").replace(tzinfo=timezone.utc).timestamp(),unit='s',tz='Europe/Budapest')
        measurement_freq = ret.getSamplingRate()  
        deltaTime = 1 / measurement_freq
        
        timestamps = pd.date_range(parsed_start_time, periods=x_acc.shape[0], freq=(str((1/measurement_freq)*1000) + "ms"))
        UFM_df = pd.DataFrame({'Timestamp': timestamps.values, 'X-Axis': x_acc, 'Y-Axis': y_acc, 'Z-Axis': z_acc, 'Magnitude': mag_acc})
        
        del x_acc, y_acc, z_acc, mag_acc, timestamps,ret
        gc.collect()
        
        fig_xyz.replace(go.Figure())
        fig_xyz.add_trace(go.Scatter(showlegend=True, name='X'), hf_x=UFM_df['Timestamp'], hf_y=UFM_df['X-Axis'], max_n_samples=5000)
        fig_xyz.add_trace(go.Scatter(showlegend=True, name='Y'), hf_x=UFM_df['Timestamp'], hf_y=UFM_df['Y-Axis'], max_n_samples=5000)
        fig_xyz.add_trace(go.Scatter(showlegend=True, name='Z'), hf_x=UFM_df['Timestamp'], hf_y=UFM_df['Z-Axis'], max_n_samples=5000)
        fig_xyz.update_layout(height=600, xaxis_title="Time", yaxis_title='Acceleration [g]', margin=dict(l=37, r=20, t=10, b=10 ), legend=dict(orientation="h"))
        
        graph_component = dcc.Graph(figure=fig_xyz, id="xyz")
        return "", graph_component

    elif file_path and string[-1] == "EDF":
        file_type = "EDF"
        signals, signal_headers, header = highlevel.read_edf(file_path, ch_nrs=selected_channels)

        x_acc = signals[0] / 1000
        y_acc = signals[1] / 1000
        z_acc = signals[2] / 1000

        global start_time_edf, patientcode, sample_rate
        start_time_edf = header['startdate']
        patientcode = header['patientcode']
        sample_rate = signal_headers[0]['sample_rate']
        deltaT = timedelta(seconds=1 / sample_rate)  
        deltaTime = deltaT.total_seconds()
        
        UFM_df = df(x_acc, y_acc, z_acc, start_time_edf, deltaT)

        del x_acc, y_acc, z_acc, signals
        gc.collect()
        
        fig_xyz.replace(go.Figure())
        fig_xyz.add_trace(go.Scatter(showlegend=True, name='X'), hf_x=UFM_df['Timestamp'], hf_y=UFM_df['X-Axis'], max_n_samples=5000)
        fig_xyz.add_trace(go.Scatter(showlegend=True, name='Y'), hf_x=UFM_df['Timestamp'], hf_y=UFM_df['Y-Axis'], max_n_samples=5000)
        fig_xyz.add_trace(go.Scatter(showlegend=True, name='Z'), hf_x=UFM_df['Timestamp'], hf_y=UFM_df['Z-Axis'], max_n_samples=5000)
        fig_xyz.update_layout(height=600, xaxis_title="Time", yaxis_title='Acceleration [g]', margin=dict(l=37, r=20, t=10, b=10 ), legend=dict(orientation="h"))
        
        graph_component = dcc.Graph(figure=fig_xyz, id="xyz")
        return "", graph_component

    return no_update, no_update

# %%
# Function to return the preprocessed DataFrame based on the method selected in the first dropdown menu
def get_dataframe(preprocessing_method):
    global UFNM_df, FMpre_df, FMpost_df, FXYZ_df, HFMpre_df
    global FX_df, FY_df, FZ_df
    global FX_df_mag, FY_df_mag, FZ_df_mag
    global UFX_df, UFY_df, UFZ_df
    global UFX_df_mag, UFY_df_mag, UFZ_df_mag
    global UFM_df, deltaTime, measurement_freq, file_type, deltaT

    if UFM_df is None:
        return None
    
    # Reset variables used in processing, only one dataframe is stored in the memory at one time
    for df_var in ['UFNM_df', 'FMpre_df', 'FMpost_df', 'FXYZ_df',
               'FX_df', 'FY_df', 'FZ_df',
               'FX_df_mag', 'FY_df_mag', 'FZ_df_mag',
               'UFX_df', 'UFY_df', 'UFZ_df',
               'UFX_df_mag', 'UFY_df_mag', 'UFZ_df_mag']:
        if df_var in globals():
            globals()[df_var] = None

    gc.collect()

    if file_type == "actigraph":
        fs = 1.0 / deltaTime
    elif file_type == "bin":
        fs = measurement_freq
    elif file_type == "EDF":
        fs = sample_rate

    
    if preprocessing_method == 'UFNM':
        UFNM_df = UFM_df.copy()
        UFNM_df['Magnitude'] = np.abs(UFNM_df['Magnitude'] - 1)


    if preprocessing_method == 'FXYZ':
        if file_type == "bin":
            FXYZ_df = pd.DataFrame({
                'Timestamp': UFM_df['Timestamp'],
                'X-Axis': butter_bandpass_filter(UFM_df['X-Axis'], fs),
                'Y-Axis': butter_bandpass_filter(UFM_df['Y-Axis'], fs),
                'Z-Axis': butter_bandpass_filter(UFM_df['Z-Axis'], fs)
            })
        else:  
            FXYZ_df = pd.DataFrame({
                'Timestamp': UFM_df['Timestamp'],
                'X-Axis': butter_bandpass_filter(UFM_df['X-Axis'], fs),
                'Y-Axis': butter_bandpass_filter(UFM_df['Y-Axis'], fs),
                'Z-Axis': butter_bandpass_filter(UFM_df['Z-Axis'], fs)
            })


    if preprocessing_method == 'FMpre':
        if file_type == "bin":
            FMpre_df = pd.DataFrame({
                'Timestamp': UFM_df['Timestamp'],
                'X-Axis': butter_bandpass_filter(UFM_df['X-Axis'], fs),
                'Y-Axis': butter_bandpass_filter(UFM_df['Y-Axis'], fs),
                'Z-Axis': butter_bandpass_filter(UFM_df['Z-Axis'], fs)
            })
            FMpre_df['Magnitude'] = np.sqrt(FMpre_df['X-Axis']**2 + FMpre_df['Y-Axis']**2 + FMpre_df['Z-Axis']**2)
        else:
            FMpre_df = pd.DataFrame({
                'Timestamp': UFM_df['Timestamp'],
                'X-Axis': butter_bandpass_filter(UFM_df['X-Axis'], fs),
                'Y-Axis': butter_bandpass_filter(UFM_df['Y-Axis'], fs),
                'Z-Axis': butter_bandpass_filter(UFM_df['Z-Axis'], fs)
            })
            FMpre_df['Magnitude'] = np.sqrt(FMpre_df['X-Axis']**2 + FMpre_df['Y-Axis']**2 + FMpre_df['Z-Axis']**2)


    if preprocessing_method == 'HFMpre':
        lowcut = 0.2
        if file_type == "bin":
            HFMpre_df = pd.DataFrame({
                'Timestamp': UFM_df['Timestamp'],
                'X-Axis': highpass_filter(UFM_df['X-Axis'], lowcut, fs),
                'Y-Axis': highpass_filter(UFM_df['Y-Axis'], lowcut, fs),
                'Z-Axis': highpass_filter(UFM_df['Z-Axis'], lowcut, fs)
            })
            HFMpre_df['Magnitude'] = np.sqrt(HFMpre_df['X-Axis']**2 + HFMpre_df['Y-Axis']**2 + HFMpre_df['Z-Axis']**2)
    
        else:
            HFMpre_df = pd.DataFrame({
                'Timestamp': UFM_df['Timestamp'],
                'X-Axis': highpass_filter(UFM_df['X-Axis'], lowcut, fs),
                'Y-Axis': highpass_filter(UFM_df['Y-Axis'], lowcut, fs),
                'Z-Axis': highpass_filter(UFM_df['Z-Axis'], lowcut, fs)
            })
            HFMpre_df['Magnitude'] = np.sqrt(HFMpre_df['X-Axis']**2 + HFMpre_df['Y-Axis']**2 + HFMpre_df['Z-Axis']**2)


    if preprocessing_method == 'FMpost':
        if file_type == "bin":
            FMpost_df = pd.DataFrame({
                'Timestamp': UFM_df['Timestamp'],
                'X-Axis': UFM_df['X-Axis'],
                'Y-Axis': UFM_df['Y-Axis'],
                'Z-Axis': UFM_df['Z-Axis']
            })
            FMpost_df['Magnitude'] = np.sqrt(FMpost_df['X-Axis']**2 + FMpost_df['Y-Axis']**2 + FMpost_df['Z-Axis']**2)
            FMpost_df['Magnitude'] = butter_bandpass_filter(FMpost_df['Magnitude'], fs)
        else:
            FMpost_df = UFM_df.copy()
            FMpost_df['Magnitude'] = np.sqrt(FMpost_df['X-Axis']**2 + FMpost_df['Y-Axis']**2 + FMpost_df['Z-Axis']**2)
            FMpost_df['Magnitude'] = butter_bandpass_filter(FMpost_df['Magnitude'], fs)


    if preprocessing_method == 'FX':
        FX_df = pd.DataFrame({
            'Timestamp': UFM_df['Timestamp'],
            'X-Axis': butter_bandpass_filter(UFM_df['X-Axis'], fs)
        })
        FX_df_mag = FX_df.rename(columns={'X-Axis': 'Magnitude'})

    if preprocessing_method == 'FY':
        FY_df = pd.DataFrame({
            'Timestamp': UFM_df['Timestamp'],
            'Y-Axis': butter_bandpass_filter(UFM_df['Y-Axis'], fs)
        })
        FY_df_mag = FY_df.rename(columns={'Y-Axis': 'Magnitude'})

    if preprocessing_method == 'FZ':
        FZ_df = pd.DataFrame({
            'Timestamp': UFM_df['Timestamp'],
            'Z-Axis': butter_bandpass_filter(UFM_df['Z-Axis'], fs)
        })
        FZ_df_mag = FZ_df.rename(columns={'Z-Axis': 'Magnitude'})


    if preprocessing_method == 'UFX':
        UFX_df = UFM_df[['Timestamp', 'X-Axis']].copy()
        UFX_df_mag = UFX_df.rename(columns={'X-Axis': 'Magnitude'})

    if preprocessing_method == 'UFY':
        UFY_df = UFM_df[['Timestamp', 'Y-Axis']].copy()
        UFY_df_mag = UFY_df.rename(columns={'Y-Axis': 'Magnitude'})

    if preprocessing_method == 'UFZ':
        UFZ_df = UFM_df[['Timestamp', 'Z-Axis']].copy()
        UFZ_df_mag = UFZ_df.rename(columns={'Z-Axis': 'Magnitude'})

    


    if preprocessing_method == 'UFNM':
        return UFNM_df
    elif preprocessing_method == 'FXYZ':
        return FXYZ_df
    elif preprocessing_method == 'FMpre':
        return FMpre_df
    elif preprocessing_method == 'FMpost':
        return FMpost_df
    elif preprocessing_method == 'HFMpre':
        return HFMpre_df
    elif preprocessing_method == 'FX':
        return FX_df
    elif preprocessing_method == 'FY':
        return FY_df
    elif preprocessing_method == 'FZ':
        return FZ_df
    elif preprocessing_method == 'UFX':
        return UFX_df
    elif preprocessing_method == 'UFY':
        return UFY_df
    elif preprocessing_method == 'UFZ':
        return UFZ_df
    else:
        return None

# %%
# Callback and function to set date pickers based on data range
@app.callback(
    Output('start-date-picker', 'min_date_allowed'),
    Output('start-date-picker', 'max_date_allowed'),
    Output('end-date-picker', 'min_date_allowed'),
    Output('end-date-picker', 'max_date_allowed'),
    Output('start-date-picker', 'date'),
    Output('end-date-picker', 'date'),
    Output('start-hour', 'value'),
    Output('start-minute', 'value'),
    Output('start-second', 'value'),
    Output('end-hour', 'value'),
    Output('end-minute', 'value'),
    Output('end-second', 'value'),
    Input('output-data-upload', 'children')
)
def update_datepickers_range(children):
    global UFM_df
    if UFM_df is not None and not UFM_df.empty:
         # Define the allowed range based on the DataFrame's timestamps
        min_date = UFM_df['Timestamp'].min()- pd.Timedelta(days=1)
        max_date = UFM_df['Timestamp'].max()- pd.Timedelta(days=1)

        start_hour = UFM_df['Timestamp'].min().hour
        start_minute = UFM_df['Timestamp'].min().minute
        start_second = UFM_df['Timestamp'].min().second
        end_hour = UFM_df['Timestamp'].max().hour
        end_minute = UFM_df['Timestamp'].max().minute
        end_second = UFM_df['Timestamp'].max().second
        if file_type == 'actigraph':
            return min_date, max_date + pd.Timedelta(days=1), min_date, max_date+ pd.Timedelta(days=1), min_date+ pd.Timedelta(days=1), max_date+ pd.Timedelta(days=1), start_hour,start_minute,start_second,end_hour,end_minute,end_second
        else:
            return min_date, max_date, min_date, max_date, min_date+ pd.Timedelta(days=1), max_date+ pd.Timedelta(days=1), start_hour,start_minute,start_second,end_hour,end_minute,end_second
    else:
        return None, None, None, None, None, None, 0,0,0,0,0,0

# %%
def generate_ufx_chart():
    fig.replace(go.Figure())
    fig.add_trace(go.Scatter(showlegend=True, name='X'), hf_x=actual_df['Timestamp'], hf_y=actual_df['X-Axis'], max_n_samples=5000)
    fig.update_layout(height=600, xaxis_title="Time", yaxis_title='Acceleration [g]', margin=dict(l=37, r=20, t=10, b=10 ), legend = dict(orientation="h"))
    return fig


def generate_ufy_chart():
    fig.replace(go.Figure())
    fig.add_trace(go.Scatter(showlegend=True, name='Y', line=dict(color='#EF553B')), hf_x=actual_df['Timestamp'], hf_y=actual_df['Y-Axis'],  max_n_samples=5000)
    fig.update_layout(height=600, xaxis_title="Time", yaxis_title='Acceleration [g]', margin=dict(l=37, r=20, t=10, b=10 ), legend = dict(orientation="h"))
    return fig


def generate_ufz_chart():
    fig.replace(go.Figure())
    fig.add_trace(go.Scatter(showlegend=True, name='Z', line=dict(color='#00CC96')), hf_x=actual_df['Timestamp'], hf_y=actual_df['Z-Axis'], max_n_samples=5000)
    fig.update_layout(height=600, xaxis_title="Time", yaxis_title='Acceleration [g]', margin=dict(l=37, r=20, t=10, b=10 ), legend = dict(orientation="h"))
    return fig


def generate_ufxyz_chart():
    fig.replace(go.Figure())
    fig.add_trace(go.Scatter(showlegend=True, name='X'), hf_x=actual_df['Timestamp'], hf_y=actual_df['X-Axis'], max_n_samples=5000)
    fig.add_trace(go.Scatter(showlegend=True, name='Y'), hf_x=actual_df['Timestamp'], hf_y=actual_df['Y-Axis'], max_n_samples=5000)
    fig.add_trace(go.Scatter(showlegend=True, name='Z'), hf_x=actual_df['Timestamp'], hf_y=actual_df['Z-Axis'], max_n_samples=5000)
    fig.update_layout(height=600, xaxis_title="Time", yaxis_title='Acceleration [g]', margin=dict(l=37, r=20, t=10, b=10 ), legend = dict(orientation="h"))
    return fig


def generate_magnitude_chart():
    fig.replace(go.Figure())
    fig.add_trace(go.Scatter(showlegend=True, name='Magnitude'), hf_x=actual_df['Timestamp'], hf_y=actual_df['Magnitude'], max_n_samples=5000)
    fig.update_layout(height=600, xaxis_title="Time", margin=dict(l=37, r=20, t=10, b=10 ), legend = dict(orientation="h"))
    return fig

# %%
def generate_TAT_chart():
    
    epoch_df = epoch(actual_df, epoch_length)
    
    global activity_df1
    if selected_preproc == "UFM":
        activity_df1 = pd.DataFrame(epoch_df['Magnitude'].groupby(epoch_df['Timestamp'], observed=False).agg(TAT, threshold=np.std(epoch_df['Magnitude'], ddof=1)+1, deltaTime=deltaTime).reset_index()).rename(columns={'Magnitude': 'Activity'})
    else:
        activity_df1 = pd.DataFrame(epoch_df['Magnitude'].groupby(epoch_df['Timestamp'], observed=False).agg(TAT, threshold=np.std(epoch_df['Magnitude'], ddof=1), deltaTime=deltaTime).reset_index()).rename(columns={'Magnitude': 'Activity'})

    fig = FigureWidgetResampler(make_subplots(specs=[[{"secondary_y": True}]]))
    generate_magnitude_chart()
    fig.add_trace(go.Scatter(showlegend=True, visible='legendonly', name='Magnitude'), hf_x=actual_df['Timestamp'], hf_y=actual_df['Magnitude'], max_n_samples=5000)
    fig.add_trace(go.Bar(showlegend=True, name='Activity', marker_line_width=0), hf_x=activity_df1['Timestamp'], hf_y=activity_df1['Activity'], max_n_samples=activity_df1['Activity'].shape[0], secondary_y=True)
    fig.update_layout(height=600, xaxis_title="Time", margin=dict(l=37, r=20, t=10, b=10 ), legend = dict(orientation="h"))
    fig.update_yaxes(title_text="Acceleration [g]", color="blue", secondary_y=False)
    fig.update_yaxes(title_text="Activity [a.u.]", color="red", secondary_y=True)

    return fig, activity_df1

def generate_MAD_chart():

    epoch_df = epoch(actual_df, epoch_length)

    global activity_df1
    activity_df1 = pd.DataFrame(epoch_df['Magnitude'].groupby(epoch_df['Timestamp'], observed=False).agg(MAD).reset_index()).rename(columns={'Magnitude': 'Activity'})
    
    fig = FigureWidgetResampler(make_subplots(specs=[[{"secondary_y": True}]]))
    generate_magnitude_chart()
    fig.add_trace(go.Scatter(showlegend=True, visible='legendonly', name='Magnitude'), hf_x=actual_df['Timestamp'], hf_y=actual_df['Magnitude'], max_n_samples=5000)
    fig.add_trace(go.Bar(showlegend=True, name='Activity', marker_line_width=0), hf_x=activity_df1['Timestamp'], hf_y=activity_df1['Activity'], max_n_samples=activity_df1['Activity'].shape[0], secondary_y=True)
    fig.update_layout(height=600, xaxis_title="Time", margin=dict(l=37, r=20, t=10, b=10), legend=dict(orientation="h"))
    fig.update_yaxes(title_text="Acceleration [g]", color="blue", secondary_y=False)
    fig.update_yaxes(title_text="Activity [a.u.]", color="red", secondary_y=True)
    
    return fig, activity_df1

def generate_ENMO_chart():

    epoch_df = epoch(actual_df, epoch_length)

    global activity_df1
    activity_df1 = pd.DataFrame(epoch_df['Magnitude'].groupby(epoch_df['Timestamp'], observed=False).agg(ENMO).reset_index()).rename(columns={'Magnitude': 'Activity'})
    
    fig = FigureWidgetResampler(make_subplots(specs=[[{"secondary_y": True}]]))
    generate_magnitude_chart()
    fig.add_trace(go.Scatter(showlegend=True, visible='legendonly', name='Magnitude'), hf_x=actual_df['Timestamp'], hf_y=actual_df['Magnitude'], max_n_samples=5000)
    fig.add_trace(go.Bar(showlegend=True, name='Activity', marker_line_width=0), hf_x=activity_df1['Timestamp'], hf_y=activity_df1['Activity'], max_n_samples=activity_df1['Activity'].shape[0], secondary_y=True)
    fig.update_layout(height=600, xaxis_title="Time", margin=dict(l=37, r=20, t=10, b=10), legend=dict(orientation="h"))
    fig.update_yaxes(title_text="Acceleration [g]", color="blue", secondary_y=False)
    fig.update_yaxes(title_text="Activity [a.u.]", color="red", secondary_y=True)
    
    return fig, activity_df1


def generate_PIM_chart():
    
    PIM_df = PIM(actual_df, epoch_length, deltaTime,selected_preproc)
    
    fig = FigureWidgetResampler(make_subplots(specs=[[{"secondary_y": True}]]))
    generate_magnitude_chart()
    fig.add_trace(go.Scatter(showlegend=True, visible='legendonly', name='Magnitude'), hf_x=actual_df['Timestamp'], hf_y=actual_df['Magnitude'], max_n_samples=5000)
    fig.add_trace(go.Bar(showlegend=True, name='Activity', marker_line_width=0), hf_x=PIM_df['Timestamp'], hf_y=PIM_df['Activity'], max_n_samples=PIM_df['Activity'].shape[0], secondary_y=True)
    fig.update_layout(height=600, xaxis_title="Time", margin=dict(l=37, r=20, t=10, b=10), legend=dict(orientation="h"))
    fig.update_yaxes(title_text="Acceleration [g]", color="blue", secondary_y=False)
    fig.update_yaxes(title_text="Activity [a.u.]", color="red", secondary_y=True)
    
    return fig, PIM_df


def generate_ZCM_chart():
    
    epoch_df = epoch(actual_df, epoch_length)

    global activity_df1
    if selected_preproc == "UFM":
        activity_df1 = pd.DataFrame(epoch_df['Magnitude'].groupby(epoch_df['Timestamp'], observed=False).agg(ZCM, threshold=np.std(epoch_df['Magnitude'], ddof=1)+1).reset_index()).rename(columns = {'Magnitude':'Activity'})
    else:
        activity_df1 = pd.DataFrame(epoch_df['Magnitude'].groupby(epoch_df['Timestamp'], observed=False).agg(ZCM, threshold=np.std(epoch_df['Magnitude'], ddof=1)).reset_index()).rename(columns = {'Magnitude':'Activity'})

    fig = FigureWidgetResampler(make_subplots(specs=[[{"secondary_y": True}]]))
    generate_magnitude_chart()
    fig.add_trace(go.Scatter(showlegend=True, visible='legendonly', name='Magnitude'), hf_x=actual_df['Timestamp'], hf_y=actual_df['Magnitude'], max_n_samples=5000)
    fig.add_trace(go.Bar(showlegend=True, name='Activity', marker_line_width=0), hf_x=activity_df1['Timestamp'], hf_y=activity_df1['Activity'], max_n_samples=activity_df1['Activity'].shape[0], secondary_y=True)
    fig.update_layout(height=600, xaxis_title="Time", margin=dict(l=37, r=20, t=10, b=10 ), legend = dict(orientation="h"))
    fig.update_yaxes(title_text="Acceleration [g]", color="blue", secondary_y=False)
    fig.update_yaxes(title_text="Activity [a.u.]", color="red", secondary_y=True)
    
    return fig, activity_df1


def generate_HFEN_chart():

    epoch_df = epoch(actual_df, epoch_length)

    global activity_df1
    activity_df1 = pd.DataFrame(epoch_df['Magnitude'].groupby(epoch_df['Timestamp'], observed=False).agg(HFEN).reset_index()).rename(columns={'Magnitude': 'Activity'})
    
    fig = FigureWidgetResampler(make_subplots(specs=[[{"secondary_y": True}]]))
    generate_magnitude_chart()
    fig.add_trace(go.Scatter(showlegend=True, visible='legendonly', name='Magnitude'), hf_x=actual_df['Timestamp'], hf_y=actual_df['Magnitude'], max_n_samples=5000)
    fig.add_trace(go.Bar(showlegend=True, name='Activity', marker_line_width=0), hf_x=activity_df1['Timestamp'], hf_y=activity_df1['Activity'], max_n_samples=activity_df1['Activity'].shape[0], secondary_y=True)
    fig.update_layout(height=600, xaxis_title="Time", margin=dict(l=37, r=20, t=10, b=10), legend=dict(orientation="h"))
    fig.update_yaxes(title_text="Acceleration [g]", color="blue", secondary_y=False)
    fig.update_yaxes(title_text="Activity [a.u.]", color="red", secondary_y=True)
    
    return fig, activity_df1

# %%
@app.callback(
    Output("graph-id", "figure"),
    Input("plot-button", "n_clicks"),
    Input("preprocessing-method", "value"),
    Input("activity-method", "value"),
    prevent_initial_call=True,
)
# Callback to update the second graph based on user inputs
def plot_graph(n_clicks, preprocessing_method, activity_method):  
    global last_df
    global actual_df
    ctx = callback_context

    if len(ctx.triggered) and "plot-button" in ctx.triggered[0]["prop_id"]:
        

        if preprocessing_method == 'UFM':
            actual_df = UFM_df
            if activity_method == 'ZCM':
                zcm_activity = generate_ZCM_chart()
                last_df = zcm_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return zcm_activity[0]
            elif activity_method == 'TAT':
                tat_activity = generate_TAT_chart()
                last_df = tat_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return tat_activity[0]
            elif activity_method == 'MAD':
                mad_activity = generate_MAD_chart()
                last_df = mad_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return mad_activity[0]
            elif activity_method == 'ENMO':
                enmo_activity = generate_ENMO_chart()
                last_df = enmo_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return enmo_activity[0]
            elif activity_method == 'PIM':
                pim_activity = generate_PIM_chart()
                last_df = pim_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return pim_activity[0]
            else:
                last_df = actual_df[['Timestamp', 'Magnitude']]
                return generate_magnitude_chart()

        elif preprocessing_method == 'UFNM':
            actual_df = get_dataframe(preprocessing_method)
            if activity_method == 'ZCM':
                zcm_activity = generate_ZCM_chart() 
                last_df = zcm_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return zcm_activity[0]
            elif activity_method == 'TAT':
                tat_activity = generate_TAT_chart()
                last_df = tat_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return tat_activity[0]
            elif activity_method == 'MAD':
                mad_activity = generate_MAD_chart()
                last_df = mad_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return mad_activity[0]
            elif activity_method == 'ENMO':
                enmo_activity = generate_ENMO_chart()
                last_df = enmo_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return enmo_activity[0]
            elif activity_method == 'PIM':
                pim_activity = generate_PIM_chart()  
                last_df = pim_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return pim_activity[0]
            else:
                last_df = actual_df[['Timestamp', 'Magnitude']]
                return generate_magnitude_chart()

        elif preprocessing_method == 'UFXYZ':  
            actual_df = UFM_df
            ufxyz_chart = generate_ufxyz_chart()
            last_df = actual_df[['Timestamp', 'X-Axis', 'Y-Axis', 'Z-Axis']]
            return ufxyz_chart
        
        elif preprocessing_method == 'UFX':
            actual_df = get_dataframe(preprocessing_method)
            if activity_method == 'MAD':
                actual_df = UFX_df_mag
                mad_activity = generate_MAD_chart()
                last_df = mad_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return mad_activity[0]
            else:
                ufx_chart = generate_ufx_chart()
                last_df = actual_df[['Timestamp', 'X-Axis']]
                last_df = last_df.rename(columns={'X-Axis': 'Magnitude'})
                return ufx_chart
        
        elif preprocessing_method == 'UFY':
            actual_df = get_dataframe(preprocessing_method)
            if activity_method == 'MAD':
                actual_df = UFY_df_mag
                mad_activity = generate_MAD_chart()
                last_df = mad_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return mad_activity[0]
            else:
                ufy_chart = generate_ufy_chart()
                last_df = actual_df[['Timestamp', 'Y-Axis']]
                last_df = last_df.rename(columns={'Y-Axis': 'Magnitude'})
                return ufy_chart
        
        elif preprocessing_method == 'UFZ':
            actual_df = get_dataframe(preprocessing_method)
            if activity_method == 'MAD':
                actual_df = UFZ_df_mag
                mad_activity = generate_MAD_chart()
                last_df = mad_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return mad_activity[0]
            else:
                ufz_chart = generate_ufz_chart()
                last_df = actual_df[['Timestamp', 'Z-Axis']]
                last_df = last_df.rename(columns={'Z-Axis': 'Magnitude'})
                return ufz_chart

        elif preprocessing_method == 'FXYZ':  
            actual_df = get_dataframe(preprocessing_method)
            fxyz_chart = generate_ufxyz_chart()
            last_df = actual_df[['Timestamp', 'X-Axis', 'Y-Axis', 'Z-Axis']]
            return fxyz_chart
            
        elif preprocessing_method == 'FX':  
            actual_df = get_dataframe(preprocessing_method)
            if activity_method == 'ZCM':
                actual_df = FX_df_mag
                zcm_activity = generate_ZCM_chart()
                last_df = zcm_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return zcm_activity[0]
            elif activity_method == 'TAT':
                actual_df = FX_df_mag
                tat_activity = generate_TAT_chart()
                last_df = tat_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return tat_activity[0]
            elif activity_method == 'MAD':
                actual_df = FX_df_mag
                mad_activity = generate_MAD_chart()
                last_df = mad_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return mad_activity[0]
            elif activity_method == 'PIM':
                actual_df = FX_df_mag
                pim_activity = generate_PIM_chart()
                last_df = pim_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return pim_activity[0]
            else:
                fx_chart = generate_ufx_chart()  
                last_df = actual_df[['Timestamp', 'X-Axis']]
                last_df = last_df.rename(columns={'X-Axis': 'Magnitude'})
                return fx_chart
            
        elif preprocessing_method == 'FY': 
            actual_df = get_dataframe(preprocessing_method)
            if activity_method == 'ZCM':
                actual_df = FY_df_mag
                zcm_activity = generate_ZCM_chart()
                last_df = zcm_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return zcm_activity[0]
            elif activity_method == 'TAT':
                actual_df = FY_df_mag
                tat_activity = generate_TAT_chart()
                last_df = tat_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return tat_activity[0]
            elif activity_method == 'MAD':
                actual_df = FY_df_mag
                mad_activity = generate_MAD_chart()
                last_df = mad_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return mad_activity[0]
            elif activity_method == 'PIM':
                actual_df = FY_df_mag
                pim_activity = generate_PIM_chart()
                last_df = pim_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return pim_activity[0]
            else:
                fy_chart = generate_ufy_chart()  
                last_df = actual_df[['Timestamp', 'Y-Axis']]
                last_df = last_df.rename(columns={'Y-Axis': 'Magnitude'})
                return fy_chart
            
        elif preprocessing_method == 'FZ':
            actual_df = get_dataframe(preprocessing_method)
            if activity_method == 'ZCM':
                actual_df = FZ_df_mag
                zcm_activity = generate_ZCM_chart()
                last_df = zcm_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return zcm_activity[0]
            elif activity_method == 'TAT':
                actual_df = FZ_df_mag
                tat_activity = generate_TAT_chart()
                last_df = tat_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return tat_activity[0]
            elif activity_method == 'MAD':
                actual_df = FZ_df_mag
                mad_activity = generate_MAD_chart()
                last_df = mad_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return mad_activity[0]
            elif activity_method == 'PIM':
                actual_df = FZ_df_mag
                pim_activity = generate_PIM_chart()
                last_df = pim_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return pim_activity[0]
            else:
                fz_chart = generate_ufz_chart()  
                last_df = actual_df[['Timestamp', 'Z-Axis']]
                last_df = last_df.rename(columns={'Z-Axis': 'Magnitude'})
                return fz_chart

        elif preprocessing_method == 'FMpre':  
            actual_df = get_dataframe(preprocessing_method)
            if activity_method == 'ZCM':
                zcm_activity = generate_ZCM_chart()
                last_df = zcm_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return zcm_activity[0]
            elif activity_method == 'TAT':
                tat_activity = generate_TAT_chart()
                last_df = tat_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return tat_activity[0]
            elif activity_method == 'MAD':
                mad_activity = generate_MAD_chart()
                last_df = mad_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return mad_activity[0]
            elif activity_method == 'PIM':
                pim_activity = generate_PIM_chart()
                last_df = pim_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return pim_activity[0]
            else:
                last_df = actual_df[['Timestamp', 'Magnitude']]
                return generate_magnitude_chart()

        elif preprocessing_method == 'FMpost':  
            actual_df = get_dataframe(preprocessing_method)
            if activity_method == 'ZCM':
                zcm_activity = generate_ZCM_chart()
                last_df = zcm_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return zcm_activity[0]
            elif activity_method == 'TAT':
                tat_activity = generate_TAT_chart()
                last_df = tat_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return tat_activity[0]
            elif activity_method == 'MAD':
                mad_activity = generate_MAD_chart()
                last_df = mad_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return mad_activity[0]
            elif activity_method == 'PIM':
                pim_activity = generate_PIM_chart()
                last_df = pim_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return pim_activity[0]
            else:
                last_df = actual_df[['Timestamp', 'Magnitude']]
                return generate_magnitude_chart()
            
        elif preprocessing_method == 'HFMpre':
            actual_df = get_dataframe(preprocessing_method)
            if activity_method == 'HFEN':
                hfen_activity = generate_HFEN_chart()
                last_df = hfen_activity[1]
                last_df = last_df.rename(columns={'Activity': 'Magnitude'})
                return hfen_activity[0]
            else:
                last_df = actual_df[['Timestamp', 'Magnitude']]
                return generate_magnitude_chart()

        else:
            return no_update
    else:
        return no_update

# %%
# Function to update the epoch length
@app.callback(
    Output('epoch-length-input', 'value'),
    Input('epoch-length-input', 'value')
)
def update_epoch_length(value):
    global epoch_length
    epoch_length = value
    return value

# Function to update the selected preprocessing method
@app.callback(
    Output('preprocessing-method', 'value'),
    [Input('preprocessing-method', 'value')]
)
def update_preprocessing_method(selected_method):
    global selected_preproc
    selected_preproc = selected_method
    return selected_method

# Function to update the bins-per-decade value for spectral analysis
@app.callback(
    Output('bins-per-decade', 'value'),  
    Input('bins-per-decade', 'value')
)
def update_bins_per_decade(bins_per_decade):
    global bins_per_decade_global
    bins_per_decade_global = bins_per_decade  
    return bins_per_decade

# Function to enable or disable buttons based on the preprocessing method
@app.callback(
    [Output('calculate-spectrum-button', 'style'),
     Output('export-csv-reduced-button', 'style')],
    [Input('preprocessing-method', 'value')]
)
def update_buttons_style(SPM):
    if SPM == 'UFXYZ' or SPM == 'FXYZ':
        disabled_style = {'backgroundColor': 'grey', 'color': 'white', 'cursor': 'not-allowed', 'pointerEvents': 'none', 'padding': '10px 20px', 'font-size': '20px', 'marginRight':'15px'}
        return disabled_style, disabled_style
    else:
        enabled_style = {'backgroundColor': '#4CAF50', 'color': '#ffffff', 'padding': '10px 20px', 'font-size': '20px', 'marginRight':'15px'}
        return enabled_style, enabled_style

# Function to toggle the activity method dropdown based on preprocessing selection
# %%
@app.callback(
    Output('activity-method', 'disabled'),
    Input('preprocessing-method', 'value')
)
def toggle_activity_dropdown(preprocessing_method):
    return preprocessing_method is None

# Function to update the activity method dropdown options based on the selected preprocessing method
@app.callback(
    Output('activity-method', 'options'),
    Input('preprocessing-method', 'value')
)
def update_second_dropdown_options(selected_value):
    if selected_value == 'UFX' or selected_value == 'UFY' or selected_value == 'UFZ' :
        updated_options = [option for option in second_dropdown_options if option['value'] in ['MAD', 'AI', '-']]
    elif selected_value == 'FX' or selected_value == 'FY' or selected_value == 'FZ':
        updated_options = [option for option in second_dropdown_options if option['value'] in ['PIM', 'ZCM', 'TAT', 'MAD', 'AI','-']]
    elif selected_value == 'UFM':
        updated_options = [option for option in second_dropdown_options if option['value'] in ['PIM', 'ZCM', 'TAT', 'MAD', 'ENMO','-']]
    elif selected_value == 'UFNM':
        updated_options = [option for option in second_dropdown_options if option['value'] in ['PIM', 'ZCM', 'TAT', 'MAD','-']]
    elif selected_value == 'FMpost' or selected_value == 'FMpre':
        updated_options = [option for option in second_dropdown_options if option['value'] in ['PIM', 'ZCM', 'TAT', 'MAD','-']]
    elif selected_value == 'HFMpre':
        updated_options = [option for option in second_dropdown_options if option['value'] in ['HFEN', '-']]
    else:
        updated_options = []
    return updated_options

# %%
# Function to export the processed data to a CSV file
@app.callback(
    Output('download-link', 'href'),
    Input('export-csv-button', 'n_clicks'),
    State('start-date-picker', 'date'),
    State('start-hour', 'value'),
    State('start-minute', 'value'),
    State('start-second', 'value'),
    State('end-date-picker', 'date'),
    State('end-hour', 'value'),
    State('end-minute', 'value'),
    State('end-second', 'value'),
    State('preprocessing-method', 'value'),
    State('activity-method', 'value'),
    prevent_initial_call=True
)
def export_to_csv(n_clicks, start_date, start_hour, start_minute, start_second, end_date, end_hour, end_minute, end_second, selected_preprocessing_method, selected_activity_metric):
    if n_clicks > 0:

        if last_df is not None:

            if start_date is None or end_date is None:
                df_data = last_df
            else:
                start_date = pd.to_datetime(start_date)
                start_datetime = start_date.replace(hour=start_hour, minute=start_minute, second=start_second)
                end_date = pd.to_datetime(end_date)
                end_datetime = end_date.replace(hour=end_hour, minute=end_minute, second=end_second)

                df_data = last_df[(last_df['Timestamp'] >= start_datetime) & (last_df['Timestamp'] < end_datetime)]


            acceleration_preprocessing_method = selected_preprocessing_method if selected_preprocessing_method else "-"
            activity_metric = selected_activity_metric if selected_activity_metric else "-"
            activity_epoch_length = epoch_length if selected_activity_metric and selected_activity_metric != "-" else "-"

            
            first_timestamp = df_data['Timestamp'].iloc[0] if not df_data.empty else "-"


            if selected_activity_metric and selected_activity_metric == "-" or selected_activity_metric == None:
                elapsed_time = deltaTime
            else:
                elapsed_time = epoch_length

            time_diff = end_datetime - start_datetime

            hours_from_days = time_diff.days * 24

            total_hours = hours_from_days + time_diff.seconds // 3600
            minutes = (time_diff.seconds % 3600) // 60
            seconds = time_diff.seconds % 60
            time_diff_formatted = f'{total_hours:02}:{minutes:02}:{seconds:02}'

            mk_part = csv_file_name[-1].split('_')[0]

            if selected_activity_metric and selected_activity_metric != "-":
                type = f"{selected_activity_metric}({acceleration_preprocessing_method})"
            else:
                type = acceleration_preprocessing_method

            def format_datetime(dt):
                if dt:
                    return dt.strftime("%d-%H-%M-%S")
                return "-"
            
            def format_datetime_year(dt):
                if dt:
                    return dt.strftime("%Y-%m-%d-%H-%M-%S")
                return "-"


            cutstart = format_datetime(start_datetime)
            cutend = format_datetime(end_datetime)

            ex_int_start = format_datetime_year(start_datetime)
            ex_int_end =  format_datetime_year(end_datetime)

            if file_type == "actigraph":
                match = re.search(r'\d{4}\.\d{2}\.\d{2}\.\d{2}', csv_file_name[-1])
                actigraph_name = match.group(0)[2:4] + match.group(0)[5:7] + match.group(0)[8:10] + match.group(0)[11:13] if match else ""

            if file_type == "EDF":
                formatted_start_time_edf = start_time_edf.strftime("%y-%m-%d")
                global patientcode
                if patientcode == '':
                    patientcode = 'X'
                else:
                    patientcode = patientcode

            if file_type == "bin":
                sampling_rate = "100"
            elif file_type == "EDF":
                sampling_rate = sample_rate
            else:
                sampling_rate = 1/deltaTime


            if file_type == "bin":
                ts_new_filename_meta = f"TS_{mk_part}_{type}_{cutstart}_{cutend}_{elapsed_time}_META.csv"
                ts_new_filename_data = f"TS_{mk_part}_{type}_{cutstart}_{cutend}_{elapsed_time}_DATA.csv"
            elif file_type == "EDF":
                ts_new_filename_meta = f"TS_{patientcode}_{formatted_start_time_edf}_{type}_{cutstart}_{cutend}_{elapsed_time}_META.csv"
                ts_new_filename_data = f"TS_{patientcode}_{formatted_start_time_edf}_{type}_{cutstart}_{cutend}_{elapsed_time}_DATA.csv"
            else:
                ts_new_filename_meta = f"TS_{actigraph_name}_{type}_{cutstart}_{cutend}_{elapsed_time}_META.csv"
                ts_new_filename_data = f"TS_{actigraph_name}_{type}_{cutstart}_{cutend}_{elapsed_time}_DATA.csv"

            data = {
                "File Name:": csv_file_name[-1],
                "Measurement Sampling Rate [Hz]:": sampling_rate,
                "Acceleration Preprocessing Method:": acceleration_preprocessing_method,
                "Activity Metric:": activity_metric,
                "Activity Epoch Length [s]:": activity_epoch_length,
                "Examination Interval Start Time:": ex_int_start,
                "Examination Interval End Time:": ex_int_end,
                "Examination Interval Length:": time_diff_formatted,
                "Timestamp of First Datapoint:": first_timestamp,
                "Elapsed Time Between Datapoints [s]:": elapsed_time 
            }

            

            df_metadata = pd.DataFrame(list(data.items()), columns=['METADATA', ''])
            df_data = df_data.drop(columns=['Timestamp'])

            # Save the data and metadata to CSV files
            csv_METADATA = os.path.join(DOWNLOAD_DIR, ts_new_filename_meta)
            csv_DATA = os.path.join(DOWNLOAD_DIR, ts_new_filename_data)
            df_data.to_csv(csv_DATA, index=False, encoding='utf-8', sep=';', header=False)
            df_metadata.to_csv(csv_METADATA, index=False, encoding='utf-8', sep=';', header=False)
            csv_base64_DATA = base64.b64encode(open(csv_DATA, 'rb').read()).decode('utf-8')
            csv_base64_METADATA = base64.b64encode(open(csv_METADATA, 'rb').read()).decode('utf-8')
            return f"data:text/csv;base64,{csv_base64_DATA}", f"data:text/csv;base64,{csv_base64_METADATA}"
    return no_update

# %%
# Function to calculate Power Spectral Density (PSD)
def calculate_psd(timestamps, magnitudes):
    df_psd = psd(timestamps, magnitudes)
    df_psd_nonzero = df_psd[df_psd['Frequency'] > 0]
    f_max = df_psd_nonzero['Frequency'].max()
    f_min = df_psd_nonzero['Frequency'].min()
    global bin_edges
    bin_edges = np.logspace(np.log10(f_min), np.log10(f_max), int(bins_per_decade_global * (np.log10(f_max) - np.log10(f_min))) + 1)
    log_bin_edges = np.log10(bin_edges)
    bin_centers = np.power(10,log_bin_edges[0:-1] + np.diff(log_bin_edges)/2)
    df_psd_reduced = pd.DataFrame({'Frequency':bin_centers, 'Magnitude':df_psd_nonzero.groupby(pd.cut(df_psd_nonzero['Frequency'], bin_edges, right=False), observed=False)['Magnitude'].mean().values}).dropna()
    return df_psd, df_psd_reduced

# %%
# Plot the spectrum graph
def plot_spectrum_graph(df_psd, df_fit, linregress_res, df_psd_reduced, fit_frequency_bound_lo=0.0001, fit_frequency_bound_hi=0.01):

    fig_spectrum.replace(go.Figure())
    fig_spectrum.add_trace(go.Scatter(showlegend=True, visible='legendonly', name='PSD', line_color="green",
                                      x=df_psd['Frequency'][1:-1], y=df_psd['Magnitude'][1:-1]),max_n_samples=100000)
    
    if linregress_res is not None:
        fig_spectrum.add_trace(go.Scatter(showlegend=True, name='Fit, exp = ' + str(np.round(linregress_res.slope, 3)), 
                                          line_width=10, line_color="blue", opacity=0.5,x=df_fit['Frequency'], y=df_fit['Magnitude']))
    
    fig_spectrum.add_trace(go.Scatter(showlegend=True, name='Averaged PSD', line_color="red", x=df_psd_reduced['Frequency'], y=df_psd_reduced['Magnitude']))
    
    for edge in bin_edges:
        fig_spectrum.add_vline(x=edge, line_color="black", opacity=0.1)

    fig_spectrum.update_layout(height=600, xaxis_title="Frequency [Hz]", yaxis_title='Magnitude [a.u.]', margin=dict(l=37, r=20, t=10, b=10), legend=dict(orientation="h"))
    fig_spectrum.update_xaxes(type="log")
    fig_spectrum.update_yaxes(type="log")
    fig_spectrum.update_xaxes(exponentformat = 'power')
    fig_spectrum.update_yaxes(exponentformat = 'power')

    return fig_spectrum

# %%
# Function to calculate and plot the spectrum
@app.callback(
    [Output('spectrum-graph', 'figure'),
     Output('slope', 'value'),
     Output('intercept', 'value'),
     Output('r2', 'value')],
    [Input('calculate-spectrum-button', 'n_clicks')],
    [State('output-data-upload', 'children'),
     State('start-date-picker', 'date'),
     State('start-hour', 'value'),
     State('start-minute', 'value'),
     State('start-second', 'value'),
     State('end-date-picker', 'date'),
     State('end-hour', 'value'),
     State('end-minute', 'value'),
     State('end-second', 'value'),
     State('fit-frequency-bound-lo', 'value'),
     State('fit-frequency-bound-hi', 'value')],
    prevent_initial_call=True
)
def calculate_and_plot_spectrum(n_clicks, children, start_date, start_hour, start_minute, start_second, 
                                end_date, end_hour, end_minute, end_second, 
                                fit_frequency_bound_lo, fit_frequency_bound_hi):
    if n_clicks > 0:
        start_date = pd.to_datetime(start_date)
        start_datetime = start_date.replace(hour=start_hour, minute=start_minute, second=start_second)
        end_date = pd.to_datetime(end_date)
        end_datetime = end_date.replace(hour=end_hour, minute=end_minute, second=end_second)

        # Filter the data based on the specified time range
        df_filtered = last_df[(last_df['Timestamp'] >= start_datetime) & (last_df['Timestamp'] <= end_datetime)]
        df_filtered = df_filtered.reset_index(drop=True)

        if not df_filtered.empty:
            global df_psd_reduced
            df_psd, df_psd_reduced = calculate_psd(df_filtered['Timestamp'], df_filtered['Magnitude'])
            
            if fit_frequency_bound_lo is not None and fit_frequency_bound_hi is not None:
                global linregress_res
                global df_fit
                # Perform linear regression on the log-log transformed PSD data
                df_fit, linregress_res = fit_line_to_psd(df_psd_reduced, fit_frequency_bound_lo, fit_frequency_bound_hi)
                spectrum_fig = plot_spectrum_graph(df_psd, df_fit, linregress_res, df_psd_reduced)

                slope = np.round(linregress_res.slope, 3)
                intercept = linregress_res.intercept
                r2 = linregress_res.rvalue ** 2

                return spectrum_fig, slope, intercept, r2
            else:
                spectrum_fig = plot_spectrum_graph(df_psd, None, None, df_psd_reduced)
                return spectrum_fig, "", "", ""

    return no_update, "", "", ""

# %%
# Function to export the averaged PSD data to a CSV
@app.callback(
    Output('download-link-reduced', 'href'),
    Input('export-csv-reduced-button', 'n_clicks'),
    State('start-date-picker', 'date'),
    State('start-hour', 'value'),
    State('start-minute', 'value'),
    State('start-second', 'value'),
    State('end-date-picker', 'date'),
    State('end-hour', 'value'),
    State('end-minute', 'value'),
    State('end-second', 'value'),
    State('preprocessing-method', 'value'),
    State('activity-method', 'value'),
    State('fit-frequency-bound-lo', 'value'),
    State('fit-frequency-bound-hi', 'value'),
    prevent_initial_call=True
)
def export_df_psd_reduced_to_csv(n_clicks,start_date, start_hour, start_minute, start_second, end_date, end_hour, end_minute, end_second, selected_preprocessing_method, selected_activity_metric, fit_frequency_bound_lo, fit_frequency_bound_hi):
    if n_clicks > 0:
        
        start_date = pd.to_datetime(start_date)
        start_datetime = start_date.replace(hour=start_hour, minute=start_minute, second=start_second)
        end_date = pd.to_datetime(end_date)
        end_datetime = end_date.replace(hour=end_hour, minute=end_minute, second=end_second)

        if selected_activity_metric and selected_activity_metric == "-" or selected_activity_metric == None:
            elapsed_time = deltaTime
        else:
            elapsed_time = epoch_length

        acceleration_preprocessing_method = selected_preprocessing_method if selected_preprocessing_method else "-"
        activity_metric = selected_activity_metric if selected_activity_metric else "-"
        activity_epoch_length = epoch_length if selected_activity_metric and selected_activity_metric != "-" else "-"

        time_diff = end_datetime - start_datetime
        hours_from_days = time_diff.days * 24
        total_hours = hours_from_days + time_diff.seconds // 3600
        minutes = (time_diff.seconds % 3600) // 60
        seconds = time_diff.seconds % 60
        time_diff_formatted = f'{total_hours:02}:{minutes:02}:{seconds:02}'

        mk_part = csv_file_name[-1].split('_')[0]

        if selected_activity_metric and selected_activity_metric != "-":
                type = f"{selected_activity_metric}({acceleration_preprocessing_method})"
        else:
                type = acceleration_preprocessing_method

        def format_datetime(dt):
                if dt:
                    return dt.strftime("%d-%H-%M-%S")
                return "-"
        
        def format_datetime_year(dt):
                if dt:
                    return dt.strftime("%Y-%m-%d-%H-%M-%S")
                return "-"

        cutstart = format_datetime(start_datetime)
        cutend = format_datetime(end_datetime)

        ex_int_start = format_datetime_year(start_datetime)
        ex_int_end =  format_datetime_year(end_datetime)

        bpd = bins_per_decade_global

        if file_type == "actigraph":
                match = re.search(r'\d{4}\.\d{2}\.\d{2}\.\d{2}', csv_file_name[-1])
                actigraph_name = match.group(0)[2:4] + match.group(0)[5:7] + match.group(0)[8:10] + match.group(0)[11:13] if match else ""

        if file_type == "EDF":
                formatted_start_time_edf = start_time_edf.strftime("%y-%m-%d")
                global patientcode
                if patientcode == '':
                    patientcode = 'X'
                else:
                    patientcode = patientcode

        if file_type == "bin":
            sampling_rate = "100"
        elif file_type == "EDF":
            sampling_rate = sample_rate
        else:
            sampling_rate = 1/deltaTime

        if file_type == 'bin':
            psd_new_filename_meta = f"PSD_{mk_part}_{type}_{cutstart}_{cutend}_{elapsed_time}_{fit_frequency_bound_lo}_{fit_frequency_bound_hi}_{bpd}_META.csv"
            psd_new_filename_data = f"PSD_{mk_part}_{type}_{cutstart}_{cutend}_{elapsed_time}_{fit_frequency_bound_lo}_{fit_frequency_bound_hi}_{bpd}_DATA.csv"
        elif file_type == "EDF":
            psd_new_filename_meta = f"PSD_{patientcode}_{formatted_start_time_edf}_{type}_{cutstart}_{cutend}_{elapsed_time}_{fit_frequency_bound_lo}_{fit_frequency_bound_hi}_{bpd}_META.csv"
            psd_new_filename_data = f"PSD_{patientcode}_{formatted_start_time_edf}_{type}_{cutstart}_{cutend}_{elapsed_time}_{fit_frequency_bound_lo}_{fit_frequency_bound_hi}_{bpd}_DATA.csv"
        else:
            psd_new_filename_meta = f"PSD_{actigraph_name}_{type}_{cutstart}_{cutend}_{elapsed_time}_{fit_frequency_bound_lo}_{fit_frequency_bound_hi}_{bpd}_META.csv"
            psd_new_filename_data = f"PSD_{actigraph_name}_{type}_{cutstart}_{cutend}_{elapsed_time}_{fit_frequency_bound_lo}_{fit_frequency_bound_hi}_{bpd}_DATA.csv"
        
        # Metadata dictionary
        data = {
            "File Name:": csv_file_name[-1],
            "Measurement Sampling Rate [Hz]:": sampling_rate,
            "Acceleration Preprocessing Method:": acceleration_preprocessing_method,
            "Activity Metric:": activity_metric,
            "Activity Epoch Length [s]:": activity_epoch_length,
            "Examination Interval Start Time:": ex_int_start,
            "Examination Interval End Time:": ex_int_end,
            "Examination Interval Length:": time_diff_formatted,
            "Bins/Decade:": bpd,
            "Fitting Interval LO Frequency [Hz]:": fit_frequency_bound_lo,
            "Fitting Interval HI Frequency [Hz]:": fit_frequency_bound_hi,
            "Slope:": np.round(linregress_res.slope, 3),
            "Intercept:": linregress_res.intercept,
            "R2:": linregress_res.rvalue**2
        }
        df = pd.DataFrame(list(data.items()), columns=['METADATA', ''])


        df_combined = pd.concat([df_psd_reduced, df_fit], axis=1)
        csv_METADATA = os.path.join(DOWNLOAD_DIR, psd_new_filename_meta)
        csv_DATA = os.path.join(DOWNLOAD_DIR, psd_new_filename_data)
        df_combined.to_csv(csv_DATA, index=False, encoding='utf-8', sep=';', header = False)
        df.to_csv(csv_METADATA, index=False, encoding='utf-8', sep=';', header=False)
        csv_base64_DATA = base64.b64encode(open(csv_DATA, 'rb').read()).decode('utf-8')
        csv_base64_METADATA = base64.b64encode(open(csv_METADATA, 'rb').read()).decode('utf-8')
        return f"data:text/csv;base64,{csv_base64_DATA}", f"data:text/csv;base64,{csv_base64_METADATA}"
    return no_update

# %%
# Register callbacks for live updates to the resampler graphs
fig_xyz.register_update_graph_callback(app=app, graph_id="xyz")
fig.register_update_graph_callback(app=app, graph_id="graph-id")
fig_spectrum.register_update_graph_callback(app=app, graph_id="spectrum-graph")

# %%
# Configuration for the app to automatically open in the browser and run on a specified port
port = 8055

def open_browser():                                                     
    webbrowser.open_new("http://localhost:{}".format(port))         

if __name__ == '__main__':
    # The next line of code is only needed when converting the .py file into .exe
    # If running in a development environment (e.g., VS Code, PyCharm), 
    # use "#" at the start of the next line to prevent opening 2 pages at the same time
    Timer(1, open_browser).start() 
    app.run(jupyter_mode="tab", port=port)

