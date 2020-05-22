# Program fft_spectrum_gui_3can_py3_01.py
# - Corrected g_scale.
# - g_scale = (Vref / ADC resolution) * (300 mv/g)
# - Error noted and corrected by Steve Ferry.
# 08/02/2018
# - Added a timeout control to a while loop.
# 12/01/2017
# Program acel_file_plot_win_01.py modified:
# - From Python 2.7 to Python 3.5.
# - Works with AVR program adxl335_3can_01.c
# 20/12/2016
# Program acel_file_plot_win_01.py
# Program fft_spectrum_gui.py modified:
# - 3 data channels (3 axes)
# 15/11/2016
# Program fft_spectrum_gui.py
# - Based on program frame_tab_plot_07.py
# - Sample acceleration data from a ADXL335 accelerometer.
# - Plot sampled data and its FFT spectrumsu
# - Save data on file and open files with saved data.
# - Serial communication with microcontroller.e
# - Serial port selection.
# - RadioButtons to select a Window function to apply.
# 13/07/2016
#hoc git hub
from datetime import datetime
import matplotlib
import os.path

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# plt.style.use('fast')
# import seaborn as sns
# sns.set(font_scale=1.5, style="whitegrid")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D, proj3d
from matplotlib import cm, dates
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.ticker as ticker

from scipy import fftpack, signal
import numpy as np
import sys

import tkinter as Tk

from tkinter import filedialog
from tkinter import messagebox
from tkinter.scrolledtext import ScrolledText
from tkinter import ttk
from tkinter import Menu
import string
import serial
import serial.tools.list_ports
import time
import math

datos_a_leer = 2048  # Amount of samples to read.16384
# datos_a_leer = 128
sample_rate = 1000  # Sampling frequency (SPS).5000
# g_scale = 3.0/512        # +- 3g. 10 bit ADC. Wrong value.
g_scale = (3.3 / 1024) * (1000 / 300)  # Right value. Thanks Steve.
# g_scale = (Vref / ADC resolution) * 1/(300 mv/g)
max_freq = sample_rate / 2  # Maximum signal frequency, X and Y axis (accelerometer).1500
max_freq_z = sample_rate / 2  # Maximum signal frequency, Z axis (accelerometer).500

g_canal_1 = []  # Global canal_1
g_canal_2 = []  # Global canal_2
g_canal_3 = []  # Global canal_3

t_timeout = 20  # Timeout time in seconds.


def rmsValue(arr):
    square = 0
    mean = 0.0
    rms = 0.0
    n = len(arr)
    # Calculate square
    for i in range(0, n):
        square += (arr[i] ** 2)

        # Calculate Mean
    mean = (square / (float)(n))
    # Calculate Root
    rms = math.sqrt(mean)
    return rms


def hamming(M):
    """Return an M + 1 point symmetric hamming window."""
    if M % 2:
        raise Exception('M must be even')
    return 0.54 - 0.46 * np.cos(2 * np.pi * np.arange(M + 1) / M)


def blackman(M):
    """Return an M + 1 point symmetric point hamming window."""
    if M % 2:
        raise Exception('M must be even')
    return (0.42 - 0.5 * np.cos(2 * np.pi * np.arange(M + 1) / M) +
            0.08 * np.cos(4 * np.pi * np.arange(M + 1) / M))


def sinc_filter(M, fc):
    """Return an M + 1 point symmetric point sinc kernel with normalised cut-off
    frequency fc 0->0.5."""
    if M % 2:
        raise Exception('M must be even')
    return np.sinc(2 * fc * (np.arange(M + 1) - M / 2))


def plot_fft(x, style='b'):
    """Convenience routine for plotting fft's of signal with normalised
    frequency axis."""
    pylab.plot(np.arange(len(x)) / float(len(x)), np.abs(np.fft.fft(x))
               , style)


def build_filter(M, fc, window=None):
    """Construct filter using the windowing method for filter parameters M
    number of taps, cut-off frequency fc and window. Window defaults to None
    i.e. a rectangular window."""
    if window is None:
        h = sinc_filter(M, fc)
    else:
        h = sinc_filter(M, fc) * window(M)
    return h / h.sum()


def scan_serial():
    """ Scans for serial ports"""
    portnames = []
    ports = (list(serial.tools.list_ports.comports()))
    for index in range(len(ports)):
        portnames.append(ports[index][0])
    return portnames


# input: string, start string, end string
# output: the string between start string and end string
def simpleParse(mainString, beginString, endString):
    """Searches for a substring between beginString and endString lay """
    posBeginString = mainString.find(beginString) + len(beginString)
    posEndString = mainString.find(endString)
    resultado = mainString[posBeginString:posEndString]
    return resultado


# input: string, preceeding tag and post tag
# output: array of number that store the sensor data
def extraer_int_tag(datos_arch, tag):
    """ Extracts data from string datos_str, delimited by <tag> y </tag>
        and convets it to integer numbers (list of integers)"""
    str_canal = ''
    beginString = '<' + tag + '>'
    endString = '</' + tag + '>'
    str_parse = simpleParse(datos_arch, beginString, endString)
    str_canal = str_parse.split(',')

    canal = []
    n = len(str_canal)
    for i in range(n):
        canal.append(float(str_canal[i]))
    return canal


def conv_str_tag(canal, tag):
    """ Convert every channel from int to str, separated by a coma
    and adds tags at the beggining and end. """
    n = len(canal)
    s_canal = '<' + tag + '>'
    for i in range(n - 1):
        s_canal = s_canal + str(canal[i]) + ','
    s_canal = s_canal + str(canal[n - 1]) + '</' + tag + '>'
    return s_canal


def grabar(canal_1, canal_2, canal_3, archivo):
    """ Saves X and Y axis data on file archivo"""
    str_canal = ''
    str_canal += conv_str_tag(canal_1, 'L1') + '\n'
    str_canal += conv_str_tag(canal_2, 'L2') + '\n'
    str_canal += conv_str_tag(canal_3, 'L3') + '\n'

    str_aux = ''
    str_aux += '<nd>' + str(datos_a_leer) + '</nd>' + '\n'
    str_aux += '<sr>' + str(sample_rate) + '<sr>' + '\n'
    # str_aux += '<gn>' + str(ganancia) + '</gn>' + '\n'

    # Write to file
    arch = open(archivo, "w")
    arch.write(str_aux)
    arch.write(str_canal)
    arch.close()


def write1channel(canal_1, archivo):
    """ Saves X and Y axis data on file archivo"""
    now = datetime.now()
    date_num = dates.date2num(now)

    str_canal = ''
    str_canal += conv_str_tag(canal_1, 'L1') + '<D>' + str(date_num) + '</D>' + '\n'

    # Write to file
    arch = open(archivo, "a")
    arch.write(str_canal)
    arch.close()


class Application():

    def __init__(self, parent):
        self.parent = parent
        self.frames()
        self.f_saved = True  # Sampled data saved
        root.protocol("WM_DELETE_WINDOW", self.on_closing)

    def on_closing(self):
        if (self.f_saved == False):
            if messagebox.askokcancel("Quit", "Sampled data not saved. Do you wanto to quit?"):
                root.destroy()
        else:
            root.destroy()

    def _quit(self):
        root.quit()
        root.destroy()
        exit()

    def frames(self):
        menuBar = Menu(root, fg="red")
        root.config(menu=menuBar)
        self.fileMenu = Menu(menuBar, tearoff=0)
        self.fileMenu1 = Menu(menuBar, tearoff=0)
        self.fileMenu2 = Menu(menuBar, tearoff=0)

        self.fileMenu.add_command(label="Import", command=self.open_file)
        self.fileMenu.add_separator()
        self.fileMenu.add_command(label="Save Data", command=self.save_file_1chanel)
        self.fileMenu.add_separator()
        self.fileMenu.add_command(label="Exit", command=self._quit)
        self.fileMenu.add_separator()

        self.fileMenu1.add_command(label="Scan Serial", command=self.scan_ports)
        self.fileMenu1.add_separator()

        self.fileMenu2.add_command(label="History Plot", command=self.open_and_read)
        self.fileMenu2.add_separator()

        menuBar.add_cascade(label="File", menu=self.fileMenu)
        menuBar.add_cascade(label="Serial", menu=self.fileMenu1)
        menuBar.add_cascade(label="Analyser", menu=self.fileMenu2)

        frame1 = Tk.Frame(root, relief='raised', borderwidth=1)
        frame2 = Tk.Frame(root, relief='raised')
        self.frame3 = Tk.Frame(root, relief='raised')
        note = ttk.Notebook(frame2)
        self.tab1 = ttk.Frame(note)
        self.tab2 = ttk.Frame(note)
        self.tab3 = ttk.Frame(note)
        self.tab4 = ttk.Frame(note)
        self.tab5 = ttk.Frame(note)

        note.add(self.tab1, text="Frequency")
        note.add(self.tab2, text="Time zone")
        note.add(self.tab3, text="History analyser")
        note.add(self.tab4, text="Monitoring")
        note.add(self.tab5, text="Balancing")
        # Positioning
        frame1.pack(side='left', fill='both')
        frame2.pack(side='left', fill='both', expand='true')

        self.frame3 = Tk.Frame(self.tab4, relief='raised')
        self.frame3.pack(side='top', fill='both', expand='true')
        self.frame4 = Tk.Frame(self.tab4, relief='raised')
        self.frame4.pack(side='bottom', fill='both', expand='true')

        self.frame5 = Tk.Frame(self.tab5, relief='raised', borderwidth=1)
        self.frame5.pack(side='left', fill='both', expand='true')
        self.frame6 = Tk.Frame(self.tab5, relief='raised', width=0.3)
        self.frame6.pack(side='right', fill='both', expand='true')

        number = Tk.StringVar()
        self.combo = ttk.Combobox(self.frame6, textvariable=number,
                                  values=["original", "1st trial mass", "2nd trial mass", "3rd trial mass"],
                                  state="readonly")
        self.combo.grid(row=0, column=0, padx=5, pady=5)
        boton_read2 = Tk.Button(self.frame6, text="READ SERIAL DATA", width=20, height=5, command=self.read_serial)
        boton_read2.grid(row=1, column=0, padx=5, pady=5)

        # boton_open = Tk.Button(frame1, text="SAVE TO FILE", width=20, height=5, command=self.save_file)
        # boton_save = Tk.Button(frame1, text="SAVE CHANNEL", width=20, height=5, command=self.save_file_1chanel)
        # boton_scan = Tk.Button(frame1, text="SCAN SERIAL PORTS", width=20, height=5, command=self.scan_ports)
        boton_read = Tk.Button(frame1, text="READ SERIAL DATA", width=20, height=5, command=self.read_serial)

        label1 = Tk.Label(frame1, text="Select Serial Port:")
        self.sel_puerto = ttk.Combobox(frame1, textvariable='', state="readonly")
        portnames = scan_serial()
        self.sel_puerto['values'] = portnames
        if (portnames != []):
            self.sel_puerto.current(0)

        self.text_message = ScrolledText(frame1, height=10, width=20, wrap=Tk.WORD)

        self.window_var = Tk.IntVar()
        self.channel_var = Tk.IntVar()
        self.filter_var = Tk.IntVar()
        self.window_var.set(1)  # Option rectangular window
        self.channel_var.set(1)
        self.filter_var.set(1)
        option = ttk.LabelFrame(frame1, text="Option")

        radio_button1 = Tk.Radiobutton(option, text="Rectangular Window",
                                       variable=self.window_var, value=1, command=self.win_sel)
        radio_button2 = Tk.Radiobutton(option, text="Hann Window",
                                       variable=self.window_var, value=2, command=self.win_sel)
        radio_button3 = Tk.Radiobutton(option, text="Flattop Window",
                                       variable=self.window_var, value=3, command=self.win_sel)

        option1 = ttk.LabelFrame(frame1, text="Choose channel before saving data")

        radio_button4 = Tk.Radiobutton(option1, text="Channel X",
                                       variable=self.channel_var, value=1)
        radio_button5 = Tk.Radiobutton(option1, text="Channel Y",
                                       variable=self.channel_var, value=2)
        radio_button6 = Tk.Radiobutton(option1, text="Channel Z",
                                       variable=self.channel_var, value=3)
        radio_button7 = Tk.Radiobutton(option1, text="All channel", value=4,
                                       variable=self.channel_var)

        option2 = ttk.LabelFrame(frame1, text="Filter-Cutoff frequency")
        radio_button8 = Tk.Radiobutton(option2, text="Lowpass",
                                       variable=self.filter_var, value=1)
        radio_button9 = Tk.Radiobutton(option2, text="Highpass",
                                       variable=self.filter_var, value=2)
        radio_button10 = Tk.Radiobutton(option2, text="Bandpass",
                                        variable=self.filter_var, value=3)

        label2 = ttk.Label(option2, text="Min")
        label3 = ttk.Label(option2, text="Max")

        self.content1 = Tk.StringVar()
        self.content2 = Tk.StringVar()
        self.content3 = Tk.StringVar()
        self.content4 = Tk.StringVar()
        text_box1 = ttk.Entry(option2, width=15, textvariable=self.content1)
        text_box2 = ttk.Entry(option2, width=15, textvariable=self.content2)
        text_box3 = ttk.Entry(option2, width=15, textvariable=self.content3)
        text_box4 = ttk.Entry(option2, width=15, textvariable=self.content4)

        button_refesh = Tk.Button(option2, text="REFRESH", width=8, height=2, command=self.change_color)
        # Grid
        # boton_save.grid(row=1, column=0, padx=5, pady=5)
        # boton_save.grid(row=2, column=0, padx=5, pady=5)
        # boton_scan.grid(row=3, column=0, padx=5, pady=5)
        label1.grid(row=0, column=0, padx=5, pady=5)
        self.sel_puerto.grid(row=1, column=0, padx=5, pady=5)
        boton_read.grid(row=2, column=0, padx=5, pady=5)
        self.text_message.grid(row=3, column=0, padx=5, pady=5)
        option.grid(row=4, column=0, sticky="W", padx=5, pady=5)
        option1.grid(row=5, column=0, sticky="W", padx=5, pady=5)
        option2.grid(row=6, column=0, sticky="W", padx=5, pady=5)

        radio_button1.grid(row=0, column=0, sticky="W")
        radio_button2.grid(row=1, column=0, sticky="W")
        radio_button3.grid(row=2, column=0, sticky="W")
        radio_button4.grid(row=0, column=0, sticky="W")
        radio_button5.grid(row=1, column=0, sticky="W")
        radio_button6.grid(row=2, column=0, sticky="W")
        radio_button7.grid(row=3, column=0, sticky="W")

        radio_button8.grid(row=0, column=0, sticky="W")
        radio_button9.grid(row=1, column=0, sticky="W")
        radio_button10.grid(row=2, column=0, sticky="W")
        label2.grid(row=3, column=0)
        label3.grid(row=4, column=0)
        text_box1.grid(row=0, column=1, sticky="W")
        text_box2.grid(row=1, column=1, sticky="W")
        text_box3.grid(row=3, column=1, sticky="W")
        text_box4.grid(row=4, column=1, sticky="W")
        button_refesh.grid(row=5, column=0, sticky="W")
        # note.grid(row = 0, column=0)
        note.pack(side='top', fill='both', padx=5, pady=5)

        # Figure 1
        fig1 = Figure(figsize=(10, 9))
        fig1.suptitle('Sampled signal - Acceleration')
        ax_11 = fig1.add_subplot(3, 1, 1)
        # ax_11.hold(False)
        ax_11.set_title("Channel X| RMS: ")
        ax_11.set_ylabel('g')
        ax_11.grid()  # Shows grid.

        ax_12 = fig1.add_subplot(3, 1, 2)
        # ax_12.hold(False)
        ax_12.set_title("Channel Y| RMS: ")
        # ax_12.set_xlabel('ms')
        ax_12.set_ylabel('g')
        ax_12.grid()  # Shows grid.

        ax_12 = fig1.add_subplot(3, 1, 3)
        # ax_12.hold(False)
        ax_12.set_title("Channel Z| RMS: ")
        ax_12.set_xlabel('ms')
        ax_12.set_ylabel('g')
        ax_12.grid()  # Shows grid.

        # Figure 2
        fig2 = Figure(figsize=(10, 9))
        fig2.suptitle('FFT spectrum')

        ax_21 = fig2.add_subplot(3, 1, 1)
        # ax_21.hold(False)
        ax_21.set_title("Channel X")
        ax_21.set_ylabel('g')
        ax_21.set_xlim(xmax=max_freq)
        ax_21.grid()

        ax_22 = fig2.add_subplot(3, 1, 2)
        # ax_22.hold(False)
        ax_22.set_title("Channel Y")
        # ax_22.set_xlabel('Hz')
        ax_22.set_xlim(xmax=max_freq)
        ax_22.set_ylabel('g')
        ax_22.grid()

        ax_23 = fig2.add_subplot(3, 1, 3)
        # ax_23.hold(False)
        ax_23.set_title("Channel Z")
        ax_23.set_xlabel('Hz')
        # ax_23.set_xlim(xmax=max_freq)
        ax_23.set_xlim(xmax=max_freq_z)
        ax_23.set_ylabel('g')
        ax_23.grid()

        # Figure 3 3D plot#
        fig3 = plt.figure()
        axis = fig3.add_subplot(1, 1, 1, projection="3d")
        fig3.subplots_adjust(left=-0.1, bottom=-0.1, right=1.04, top=1.07)
        axis.view_init(azim=310, elev=15)
        # 290,15
        axis.set_xlabel("X")
        axis.set_ylabel('Y')
        axis.set_zlabel('Z')
        # axis.set_title('HISTORY')

        # the number of file store in the path
        # img_folder_path = 'C:/Users/Admin/Desktop/python/history data/'
        # dirListing = os.listdir(img_folder_path)

        gridsize = (3, 2)
        fig4 = plt.figure(figsize=(12, 10))
        ax1 = plt.subplot2grid(gridsize, (1, 0), colspan=2, rowspan=2)
        # ax2 = plt.subplot2grid(gridsize, (0, 0))
        ax3 = plt.subplot2grid(gridsize, (0, 1))
        fig4.subplots_adjust(bottom=0.1, right=0.9, top=0.95, hspace=0.3)

        # rms_arr2 = 1.4 * np.ones(3)
        # rms_arr3 = 2.8 * np.ones(3)
        # rms_arr4 = 4.5 * np.ones(3)
        # time_table=[1,2,3]
        # ax1.plot(time_table, rms_arr2, color="green")
        # ax1.plot(time_table, rms_arr3, color="yellow")
        # ax1.plot(time_table, rms_arr4, color="red")
        ax1.set_title("RMS BELONG TIME")
        ax1.set_xlabel('time')
        ax1.set_ylabel('Standard Velocity RMS value ISO-10816 (mm/s)')
        labels = 'MISALIGMENT', 'BELT', 'UNBALANCE', 'LOOSE', 'BENT SHAFT'
        sizes = [20, 20, 20, 20, 20]
        explode = (0.1, 0, 0, 0, 0)  # only "explode" the 2nd slice (i.e. 'Hogs')
        ax3.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%',
                shadow=True, startangle=90)
        ax3.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
        ax3.set_title('Fault type')
        xs = [0, 3.14]
        org = [2.6, 6.5, 1.9]
        trial = org[1:3]
        vt = np.sqrt((org[1] ** 2 + org[2] ** 2 - 2 * (org[0] ** 2)) / 2)
        amp = [vt, vt]
        theta = np.arccos((org[1] ** 2 - org[2] ** 2) / (4 * vt * org[0]))
        mcom = org[0] * 10 / vt
        # print(trial)
        # print("vt=", vt)
        # print('theta=', theta)
        fig5 = plt.figure(figsize=(11, 11))
        ax_51 = plt.subplot(2, 1, 1, projection='polar', autoscale_on=True)
        for x1, y1 in zip(xs, amp):
            plt.polar(x1, y1, "*")
            plt.text(x1, y1, '(%0.2f, %0.1f)' % (x1, y1))

        plt.polar(theta, org[0], "*")
        plt.text(theta, org[0], '(%0.2f, %0.1f)' % (theta, org[0]))

        plt.polar(theta + np.pi, org[0], "*")
        plt.text(theta + np.pi, org[0], '(comp mass %0.2f, %0.1f, %0.1f)' % (theta + np.pi, org[0], mcom))
        plt.plot([0, theta + np.pi], [0, org[0]], '-y')
        plt.plot([0, 0], [0, vt], '-g')
        plt.plot([0, np.pi], [0, vt], '-b')
        plt.plot([0, theta], [0, org[0]], '-r')
        # plt.subplots_adjust(top=0.95, bottom=0.0, left=0.001, right=1,autoscale_on=True)

        self.canvas5 = FigureCanvasTkAgg(fig5, master=self.frame5)
        self.canvas5.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        self.toolbar5 = NavigationToolbar2Tk(self.canvas5, self.frame5)
        self.toolbar5.update()
        self.canvas5._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        self.canvas4 = FigureCanvasTkAgg(fig4, master=self.frame3)
        self.canvas4.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        self.toolbar4 = NavigationToolbar2Tk(self.canvas4, self.frame3)
        self.toolbar4.update()
        self.canvas4._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        self.canvas3 = FigureCanvasTkAgg(fig3, master=self.tab3)
        self.canvas3.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        self.toolbar3 = NavigationToolbar2Tk(self.canvas3, self.tab3)
        self.toolbar3.update()
        self.canvas3._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        self.canvas2 = FigureCanvasTkAgg(fig1, master=self.tab2)
        self.canvas2.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        self.toolbar2 = NavigationToolbar2Tk(self.canvas2, self.tab2)
        self.toolbar2.update()
        self.canvas2._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        self.canvas1 = FigureCanvasTkAgg(fig2, master=self.tab1)
        self.canvas1.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        self.toolbar1 = NavigationToolbar2Tk(self.canvas1, self.tab1)
        self.toolbar1.update()
        self.canvas1._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

    def read_serial(self):
        puerto = self.sel_puerto.get()
        print(puerto)
        message_string = "Port: {0} \n".format(puerto)
        self.show_message(self.text_message, message_string)

        estado_serial = False
        try:
            serial_avr = serial.Serial(port=puerto, baudrate=500000,
                                       bytesize=serial.EIGHTBITS, parity=serial.PARITY_NONE,
                                       stopbits=serial.STOPBITS_ONE, timeout=0)

            time.sleep(2)  # waiting the initialization...
            print("Initializing")
            message_string = "Initializing... \n"
            self.show_message(self.text_message, message_string)

            if (serial_avr.isOpen() == True):
                estado_serial = True
            else:
                estado_serial = False
        except (serial.SerialException, ValueError) as ex:
            # print "CanÂ´t open serial port: " + str(ex)
            messagebox.showerror("Result", "Can't open serial port: " + str(ex))

        if (estado_serial == True):
            global g_canal_1, g_canal_2, g_canal_3, datos_a_leer
            canal_1 = []
            canal_2 = []
            canal_3 = []
            buffer = []
            paquete = []
            valores = []
            serial_avr.flushInput()  # xoa bo nho dem
            serial_avr.flushOutput()

            valores_decod1 = []
            valores_decod2 = []
            valores_decod3 = []
            defalt_val1 = [1, 0, 0, 0, 0, 0]
            defalt_val2 = [1, 0, 0, 0, 0, 0]
            defalt_val3 = [1, 0, 0, 0, 0, 0]

            conta_datos_rx = 0;  # Received samples counter.

            print("Sending INI")
            message_string = "Sending INI \n"
            self.show_message(self.text_message, message_string)

            serial_avr.write(b'INI')  # Start data sampling command.
            # serial_avr.write(chr(0x22))    #CRC 'INI'. Not used.
            serial_avr.write(b"\x7E")  # End of packet.
            time.sleep(2)
            global t_timeout
            timeout_state = False
            t0 = time.time()  # Start loop time stamp.
            while ((conta_datos_rx < datos_a_leer) and (timeout_state == False)):

                if serial_avr.inWaiting():
                    # print(serial_avr.inWaiting())
                    lectura = serial_avr.read(serial_avr.inWaiting())  # it return the how many byte containing in the buffer
                    buffer += lectura
                    valores = []
                if len(buffer) > 11:
                    try:
                        i = buffer.index(0x7E)  # tim ki tu 0x7E tra ve vi tri cua no
                    except (ValueError):
                        i = -1
                    # print("Buffer: {0}".format(buffer))
                    if i >= 0:
                        paquete = buffer[:i]  # lay i phan tu dau
                        buffer = buffer[i + 1:]  # lay het tu phan tu thu i+1
                        # print("Paquete: {0}".format(paquete))
                        valores = [i for i in paquete]
                        # print(valores)
                        # print(len(valores))
                        paquete = ''
                        if (len(valores) == 18):
                            x = 0
                            count = 0
                            # try:
                            while x < 15:
                                if (((valores[x] == 0x2D) or (valores[x] == 0x2B)) and (count == 0)):
                                    if (valores[x] == 0x2D):
                                        valores_decod1.append(-1)
                                    else:
                                        valores_decod1.append(1)
                                    valores_decod1.append(valores[x + 1] - 48)
                                    valores_decod1.append(valores[x + 2] - 48)
                                    valores_decod1.append(valores[x + 3] - 48)
                                    valores_decod1.append(valores[x + 4] - 48)
                                    valores_decod1.append(valores[x + 5] - 48)
                                    x = x + 1
                                    count += 1
                                    # print(x)
                                if (((valores[x] == 0x2D) or (valores[x] == 0x2B)) and (count == 1)):
                                    if (valores[x] == 0x2D):
                                        valores_decod2.append(-1)
                                    else:
                                        valores_decod2.append(1)
                                    valores_decod2.append(valores[x + 1] - 48)
                                    valores_decod2.append(valores[x + 2] - 48)
                                    valores_decod2.append(valores[x + 3] - 48)
                                    valores_decod2.append(valores[x + 4] - 48)
                                    valores_decod2.append(valores[x + 5] - 48)

                                    x = x + 1
                                    count += 1
                                    # print(x)
                                if (((valores[x] == 0x2D) or (valores[x] == 0x2B)) and (count == 2)):
                                    # print(x)
                                    if (valores[x] == 0x2D):
                                        valores_decod3.append(-1)
                                    else:
                                        valores_decod3.append(1)
                                    valores_decod3.append(valores[x + 1] - 48)
                                    valores_decod3.append(valores[x + 2] - 48)
                                    valores_decod3.append(valores[x + 3] - 48)
                                    valores_decod3.append(valores[x + 4] - 48)
                                    valores_decod3.append(valores[x + 5] - 48)
                                    x = x + 1
                                    count = 0
                                    # print(x)
                                else:
                                    # valores_decod.append(valores[x])
                                    x = x + 1
                            # except IndexError:
                            #     print(x)
                            #     print(valores)
                            # print("val1:", valores_decod1)
                            # print("val2:", valores_decod2)
                            # print("val3:", valores_decod3)
                        else:
                            valores_decod1.append(defalt_val1[0])
                            valores_decod1.append(defalt_val1[1])
                            valores_decod1.append(defalt_val1[2])
                            valores_decod1.append(defalt_val1[3])
                            valores_decod1.append(defalt_val1[4])
                            valores_decod1.append(defalt_val1[5])

                            valores_decod2.append(defalt_val2[0])
                            valores_decod2.append(defalt_val2[1])
                            valores_decod2.append(defalt_val2[2])
                            valores_decod2.append(defalt_val2[3])
                            valores_decod2.append(defalt_val2[4])
                            valores_decod2.append(defalt_val2[5])

                            valores_decod3.append(defalt_val3[0])
                            valores_decod3.append(defalt_val3[1])
                            valores_decod3.append(defalt_val3[2])
                            valores_decod3.append(defalt_val3[3])
                            valores_decod3.append(defalt_val3[4])
                            valores_decod3.append(defalt_val3[5])
                            print("chieu dai doan buffer nhan duoc:",len(valores))
                            print(conta_datos_rx)
                        conta_datos_rx += 1
                        # print(conta_datos_rx)
                        canal1 = valores_decod1[0] * (
                                    (valores_decod1[1] * 10000) + valores_decod1[2] * 1000 + valores_decod1[3] * 100 +
                                    valores_decod1[4] * 10 + valores_decod1[5]) / 1000
                        canal2 = valores_decod2[0] * (
                                    (valores_decod2[1] * 10000) + valores_decod2[2] * 1000 + valores_decod2[3] * 100 +
                                    valores_decod2[4] * 10 + valores_decod2[5]) / 1000
                        canal3 = valores_decod3[0] * (
                                    (valores_decod3[1] * 10000) + valores_decod3[2] * 1000 + valores_decod3[3] * 100 +
                                    valores_decod3[4] * 10 + valores_decod3[5]) / 1000

                        canal_1.append(canal1)
                        canal_2.append(canal2)
                        canal_3.append(canal3)

                        print("Canal 1: %s    Canal2: %s  " % (canal1, canal2))
                        defalt_val1 = valores_decod1
                        defalt_val2 = valores_decod2
                        defalt_val3 = valores_decod3
                        valores = []
                        valores_decod1 = []
                        valores_decod2 = []
                        valores_decod3 = []

                # conta_datos_rx += 1;
                # print("conta_datos_rx =  %s" %conta_datos_rx)

                # Check if t_timeout seconds have elapsed since time stamp t0
                if ((time.time() - t0) > t_timeout):
                    timeout_state = True
                    # print("Serial port timeout")

            if (timeout_state == False):
                print(time.time() - t0)
                print("Sending PAR")
                self.text_message.config(state=Tk.NORMAL)  # Enable to modify
                self.text_message.insert(Tk.END, "Sending PAR \n")
                self.text_message.config(state=Tk.DISABLED)  # Disable - Read only
                root.update_idletasks()  # Needed to make message visible

                serial_avr.write(b'PAR')  # Stop data sampling.
                serial_avr.write(b"\x7E")  # End of packet.

                serial_avr.close()  # Close serial port.

                print("Amount of samples channel 1: %s" % len(canal_1))
                print("Amount of samples channel 2: %s" % len(canal_2))
                print("Amount of samples channel 3: %s" % len(canal_3))
                message_string = "Amount of samples channel 1: {0} \n".format(len(canal_1))
                message_string += "Amount of samples channel 2: {0} \n".format(len(canal_2))
                message_string += "Amount of samples channel 3: {0} \n".format(len(canal_3))
                self.show_message(self.text_message, message_string)

                # Keep a copy of the original values
                g_canal_1 = canal_1[:]  # Copy list by value not by reference
                g_canal_2 = canal_2[:]
                g_canal_3 = canal_3[:]

                self.f_saved = False  # Sampled data not saved

                self.window_var.set(1)  # Option rectangular window
                self.plot(self.tab1, self.tab2, canal_1, canal_2, canal_3, win_var=1)
            else:
                serial_avr.write(b'PAR')  # Stop data sampling.
                serial_avr.write(b"\x7E")  # End of packet.
                serial_avr.close()  # Close serial port.
                print("Serial port timeout")
                message_string = ("Serial port timeout \n")
                self.show_message(self.text_message, message_string)

    def show_message(self, text_message, message_string):
        """Shows messages on a scrollable textbox"""
        text_message.config(state=Tk.NORMAL)  # Enable to modify
        text_message.insert(Tk.END, message_string)
        text_message.config(state=Tk.DISABLED)  # Disable - Read only
        text_message.see("end")  # Show the "end" of text
        root.update_idletasks()  # Needed to make message visible

    def scan_ports(self):
        portnames = []
        portnames = scan_serial()
        self.fileMenu1.entryconfigure(0, label="Scan Serial" + "(" + str(portnames) + ")")
        self.sel_puerto['values'] = portnames
        if (portnames != []):
            self.sel_puerto.current(0)

    def plot(self, tab1, tab2, canal_1, canal_2, canal_3, win_var=1):
        num_datos = len(canal_1)
        X = range(0, num_datos, 1)

        # Scale the signal in g's
        # for indice in X:
        #     canal_1[indice] *= g_scale
        #     canal_2[indice] *= g_scale
        #     canal_3[indice] *= g_scale

        # Calculates medium value for each channel.
        vdc_canal_1 = 0
        vdc_canal_2 = 0
        vdc_canal_3 = 0
        for indice in X:
            vdc_canal_1 += canal_1[indice]
            vdc_canal_2 += canal_2[indice]
            vdc_canal_3 += canal_3[indice]

        vdc_canal_1 = vdc_canal_1 / num_datos
        vdc_canal_2 = vdc_canal_2 / num_datos
        vdc_canal_3 = vdc_canal_3 / num_datos
        # print("Vdc Channel 1: {0}, Vdc Channel 2: {1}".format(vdc_canal_1, vdc_canal_2))
        # dua ve gia tri trung binh 0
        # Substract DC offset
        for indice in X:
            canal_1[indice] -= 0
            canal_2[indice] -= 0
            canal_3[indice] -= 0

        # ----------------- Plotting ----------
        X1 = np.linspace(0, int(1.85 * num_datos), num=num_datos)  # X axis, 5000 sps, 1/5 ms.

        # Figure 1. Sampled signals.
        # Channel X
        rms_valX = rmsValue(canal_1)
        rms_valY = rmsValue(canal_2)
        rms_valZ = rmsValue(canal_3)
        titX = "Channel X" + "|RMS: " + str(rms_valX)[:5]
        titY = "Channel Y" + "|RMS: " + str(rms_valY)[:5]
        titZ = "Channel Z" + "|RMS: " + str(rms_valZ)[:5]
        ax_11, ax_12, ax_13 = self.canvas2.figure.get_axes()
        ax_11.clear()
        ax_11.plot(X1, canal_1, color='blue', linewidth=0.25)
        ax_11.set_title(titX)
        ax_11.set_ylabel('g')
        ax_11.grid()  # Shows grid.

        # Channel Y
        ax_12.clear()
        ax_12.plot(X1, canal_2, color='blue', linewidth=0.25)
        ax_12.set_title(titY)
        # ax_12.set_xlabel('ms')
        ax_12.set_ylabel('g')
        ax_12.grid()  # Shows grid.

        # Channel Z
        ax_13.clear()
        ax_13.plot(X1, canal_3, color='blue', linewidth=0.25)
        ax_13.set_title(titZ)
        ax_13.set_xlabel('ms')
        ax_13.set_ylabel('g')
        ax_13.grid()  # Shows grid.

        # Figure 2. FFT from signals.
        # Channel X
        canal_fft = []
        canal_fft = canal_1

        N = len(canal_fft)  # length of the signal

        # Window function
        if (win_var == 2):
            w = signal.hann(N, sym=False)  # Hann (Hanning) window
            # print("wh=",len(w))
        elif (win_var == 3):
            w = signal.flattop(N, sym=False)  # Flattop window
            # print("wf=",len(w))
        else:
            w = 1  # Rectangular window

        T = 1.0 / sample_rate
        y = canal_fft
        yf = fftpack.fft(y * w) * (2 / N)
        yf = yf[:(int(N / 2))]
        # kenh1=yf
        xf = np.linspace(0.0, 1.0 / (2.0 * T), int(N / 2))
        # print(xf[:10])
        ax_21, ax_22, ax_23 = self.canvas1.figure.get_axes()
        ax_21.clear()
        # ax_21.plot(xf, 2.0/N * np.abs(yf[:N/2]))
        ax_21.plot(xf, np.abs(yf), color='blue', linewidth=0.5)
        ax_21.grid()
        ax_21.set_title("Channel X")
        ax_21.set_ylabel('g')
        ax_21.set_xlim(xmax=max_freq)

        # Channel Y
        canal_fft = []
        canal_fft = canal_2

        N = len(canal_fft)  # length of the signal
        T = 1.0 / sample_rate
        y = canal_fft
        yf = fftpack.fft(y * w) * (2 / N)
        yf = yf[:int(N / 2)]
        # kenh2=yf
        xf = np.linspace(0.0, 1.0 / (2.0 * T), int(N / 2))

        ax_22.clear()
        # ax_22.plot(xf, 2.0/N * np.abs(yf[:N/2]))
        ax_22.plot(xf, np.abs(yf), color='blue', linewidth=0.5)
        ax_22.grid()
        ax_22.set_title("Channel Y")
        # ax_22.set_xlabel('Hz')
        ax_22.set_xlim(xmax=max_freq)
        ax_22.set_ylabel('g')

        # Channel Z
        canal_fft = []
        canal_fft = canal_3

        N = len(canal_fft)  # length of the signal
        T = 1.0 / sample_rate
        y = canal_fft
        yf = fftpack.fft(y * w) * (2 / N)
        yf = yf[:int(N / 2)]
        # kenh3=yf
        xf = np.linspace(0.0, 1.0 / (2.0 * T), int(N / 2))
        ax_23.clear()
        # ax_23.plot(xf, 2.0/N * np.abs(yf[:N/2]))
        ax_23.plot(xf, np.abs(yf), color='blue', linewidth=0.5)
        ax_23.grid()
        ax_23.set_title("Channel Z")
        ax_23.set_xlabel('Hz')
        # ax_23.set_xlim(xmax=max_freq)
        ax_23.set_xlim(xmax=max_freq_z)
        ax_23.set_ylabel('g')
        self.canvas1.draw()
        self.canvas2.draw()
        # self.canvas3.draw()

    def plot_all_history(self, tab3, tab4, arr, d_arr, win_var=1):
        h_num = len(arr)
        # h_arr=[]
        N = len(arr[0])  # length of the signal
        T = 1.0 / sample_rate
        # Window function
        rms_arr = []
        rms_arr2 = 1.4 * np.ones(h_num)
        rms_arr3 = 2.8 * np.ones(h_num)
        rms_arr4 = 4.5 * np.ones(h_num)
        date_arr = []
        for i in range(0, h_num):
            rms_arr.append(rmsValue(arr[i]))

            date_arr.append(d_arr[i][0])
        if (win_var == 2):
            w = signal.hann(N, sym=False)  # Hann (Hanning) window
            # print("wh=",len(w))
        elif (win_var == 3):
            w = signal.flattop(N, sym=False)  # Flattop window
            # print("wf=",len(w))
        else:
            w = 1  # Rectangular window

        ax_31, = self.canvas3.figure.get_axes()
        ax_31.clear()
        # ax_31.set_title("HISTORY ANALYSER")
        ax_31.set_xlabel('Frequence (Hz)')
        ax_31.set_ylabel('Time')
        ax_31.set_zlabel('Ampli(g)')
        ax_31.grid(row=0, column=0, sticky="W", padx=1, pady=1)
        xf = np.linspace(0.0, 1.0 / (2.0 * T), int(N / 2))
        zf = []
        hfmt = dates.DateFormatter('%Y-%m-%d')
        # for h_count in range(0,h_num):
        #     temp_zf=d_arr[h_count][0]*np.ones(int(N/2))
        #     zf.append(date_num)
        #     print(date_num)
        for i in range(0, h_num):
            y = arr[i]
            h_sum = 0
            for h_count in range(0, N):
                h_sum += y[h_count]
            h_sum = h_sum / N
            for h_count in range(0, N):
                y[h_count] -= h_sum
            yf = fftpack.fft(y * w) * (2 / N)
            yf = yf[:int(N / 2)]
            try:
                zz = d_arr[i][0] * np.ones(len(xf))
                # print(zz)
                # print(len(xf))
            except:
                pass
            # zf = 15*np.ones(len(xf))
            ax_31.plot3D(xf, zz, np.abs(yf), linewidth=0.5)
            ax_31.yaxis.set_major_locator(dates.AutoDateLocator())
            ax_31.yaxis.set_major_formatter(hfmt)
        ax_41, ax_43 = self.canvas4.figure.get_axes()
        ax_41.clear()
        ax_43.clear()
        ax_41.plot(date_arr, rms_arr, color="blue", label="RMS Value")
        ax_41.plot(date_arr, rms_arr2, color="green", label="GOOD LEVEL")
        ax_41.plot(date_arr, rms_arr3, color="yellow", label="WARNING LEVEL")
        ax_41.plot(date_arr, rms_arr4, color="red", label="DAMAGE LEVEL")
        ax_41.set_title("RMS BELONG TIME")
        ax_41.set_xlabel('time')
        ax_41.set_ylabel('Standard Velocity RMS value ISO-10816 (mm/s)')
        ax_41.xaxis.set_major_locator(dates.AutoDateLocator())
        ax_41.xaxis.set_major_formatter(hfmt)
        ax_41.legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0.)
        labels = 'MISALIGMENT', 'BELT', 'UNBALANCE', 'LOOSE', 'BENT SHAFT'
        sizes = [67, 20, 10, 5, 5]
        explode = (0.1, 0.1, 0.1, 0.1, 0.1)  # only "explode" the 2nd slice (i.e. 'Hogs')
        ax_43.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%',
                  shadow=True, startangle=90)
        ax_43.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
        ax_43.set_title('Fault type')
        self.canvas3.draw()
        self.canvas4.draw()

    def win_sel(self):
        """Window selection. Every time a window is selected,
        the FFT spectrum is calculated, applying the selected window function"""
        global g_canal_1, g_canal_2, g_canal_3
        canal_1 = g_canal_1[:]  # Copy list by value not by reference
        canal_2 = g_canal_2[:]
        canal_3 = g_canal_3[:]
        win_var = self.window_var.get()
        if (len(canal_1) != 0):  # Apply only if data available
            self.plot(self.tab1, self.tab2, canal_1, canal_2, canal_3, win_var)

    def open_file(self):
        """Opens dialog to select a file, reads data from file and plots the data"""
        ftypes = [('Text files', '*.txt'), ('All files', '*')]
        dlg = filedialog.Open(root, filetypes=ftypes)
        fl = dlg.show()
        if fl != '':
            # Open file for reading
            arch = open(fl, "r")
            datos_arch = arch.read()
            # Searches for every channel, delimited by L1, L2 and L3 tags.
            canal_1 = extraer_int_tag(datos_arch, 'L1')
            canal_2 = extraer_int_tag(datos_arch, 'L2')
            canal_3 = extraer_int_tag(datos_arch, 'L3')

            print("Amount of samples in channel 1: %s" % len(canal_1))
            print("Amount of samples on channel 2: %s" % len(canal_2))
            print("Amount of samples on channel 3: %s" % len(canal_3))

            message_string = "Amount of samples channel 1: {0} \n".format(len(canal_1))
            message_string += "Amount of samples channel 2: {0} \n".format(len(canal_2))
            message_string += "Amount of samples channel 3: {0} \n".format(len(canal_3))
            self.show_message(self.text_message, message_string)

            global g_canal_1, g_canal_2, g_canal_3
            # Keep a copy of the original values
            g_canal_1 = canal_1[:]  # Copy list by value not by reference
            g_canal_2 = canal_2[:]
            g_canal_3 = canal_3[:]

            self.window_var.set(1)  # Option rectangular window
            self.plot(self.tab1, self.tab2, canal_1, canal_2, canal_3, win_var=1)

    def open_and_read(self):
        pass
        ftypes = [('Text files', '*.txt'), ('All files', '*')]
        dlg = filedialog.Open(root, filetypes=ftypes)
        fl1 = dlg.show()
        if fl1 != '':
            file1 = open(fl1, "r")
            data = file1.readlines()
            # print(data[3])
            array1 = []
            date_arr = []
            for line in range(2, len(data)):
                temp_arr = extraer_int_tag(data[line], 'L1')
                temp_date = extraer_int_tag(data[line], 'D')
                array1.append(temp_arr)
                date_arr.append(temp_date)
            self.plot_all_history(self.tab3, self.tab4, array1, date_arr, win_var=1)

    def save_file(self):
        ftypes = [('Text files', '*.txt'), ('All files', '*')]
        dlg = filedialog.SaveAs(root, filetypes=ftypes)
        fl = dlg.show()
        if fl != '':
            global g_canal_1, g_canal_2, g_canal_3
            if (len(g_canal_1) > 0):
                grabar(g_canal_1, g_canal_2, g_canal_3, fl)
                self.f_saved = True  # Sampled data saved
            else:
                print("No sampled data to save")
                message_string = "No sampled data to save\n"
                self.show_message(self.text_message, message_string)

    def save_file_1chanel(self):
        ftypes = [('Text files', '*.txt'), ('All files', '*')]
        dlg = filedialog.SaveAs(root, filetypes=ftypes)
        fl = dlg.show()
        if fl != '':
            # fl='C:/Users/Admin/Desktop/python/history data/company2/otani_analys_history_chanelX.txt'
            global g_canal_1, g_canal_2, g_canal_3
            channelVal = self.channel_var.get()
            if (len(g_canal_1) > 0):
                if (channelVal == 1):
                    write1channel(g_canal_1, fl)
                    self.f_saved = True  # Sampled data saved
                elif (channelVal == 2):
                    write1channel(g_canal_2, fl)
                    self.f_saved = True  # Sampled data saved
                elif (channelVal == 3):
                    write1channel(g_canal_3, fl)
                    self.f_saved = True
                else:
                    grabar(g_canal_1, g_canal_2, g_canal_3, fl)
                    self.f_saved = True

            else:
                print("No sampled data to save")
                message_string = "No sampled data to save\n"
                self.show_message(self.text_message, message_string)

    def change_color(self):

        global g_canal_1, g_canal_2, g_canal_3
        canal_1 = g_canal_1[:]  # Copy list by value not by reference
        canal_2 = g_canal_2[:]
        canal_3 = g_canal_3[:]

        filter_val = self.filter_var.get()
        M = int(3.21 / (10 / sample_rate))
        if (M % 2 == 1):
            M += 1
        shift = np.cos(2 * np.pi * 0.5 * np.arange(M + 1))

        if (filter_val == 1):
            temp_content = int(self.content1.get())
            fc = temp_content / sample_rate
            ham_lp = build_filter(M, fc, window=hamming)
            f_ham1 = np.convolve(canal_1, ham_lp)
            f_ham2 = np.convolve(canal_2, ham_lp)
            f_ham3 = np.convolve(canal_3, ham_lp)
            print(len(f_ham1[M:16384]))
            if (len(f_ham1) != 0):  # Apply only if data available
                self.plot(self.tab1, self.tab2, f_ham1[M:16384], f_ham2[M:16384], f_ham3[M:16384], win_var=1)
        elif (filter_val == 2):
            temp_content = int(self.content2.get())
            fc = 0.5 - (temp_content / sample_rate)
            ham_lp = build_filter(M, fc, window=hamming)
            ham_hp = ham_lp * shift
            f_ham1 = np.convolve(canal_1, ham_hp)
            f_ham2 = np.convolve(canal_2, ham_hp)
            f_ham3 = np.convolve(canal_3, ham_hp)
            if (len(f_ham1) != 0):  # Apply only if data available
                self.plot(self.tab1, self.tab2, f_ham1[M:16384], f_ham2[M:16384], f_ham3[M:16384], win_var=1)
        else:
            temp_content = int(self.content3.get())
            temp_content1 = int(self.content4.get())
            fc2 = 0.5 - (temp_content / sample_rate)
            fc1 = (temp_content1 / sample_rate)
            ham_lp = build_filter(M, fc1, window=hamming)
            ham_lp1 = build_filter(M, fc2, window=hamming)
            ham_hp = ham_lp1 * shift
            ham_bp = np.convolve(ham_lp, ham_hp)
            f_ham1 = np.convolve(canal_1, ham_bp)
            f_ham2 = np.convolve(canal_2, ham_bp)
            f_ham3 = np.convolve(canal_3, ham_bp)
            print(len(f_ham1))
            if (len(f_ham1) != 0):  # Apply only if data available
                self.plot(self.tab1, self.tab2, f_ham1[2 * M:16384], f_ham2[2 * M:16384], f_ham3[2 * M:16384],
                          win_var=1)


if __name__ == '__main__':
    root = Tk.Tk()
    # root.wm_attributes('-fullscreen','true')
    # root.attributes('-fullscreen', True)
    # root.call("wm", "attributes", ".", "-fullscreen", "true")  # Fullscreen mode
    root.title('OTANI-UP FFT ANALYSER')
    app = Application(root)
    root.mainloop()
