# DMREF_Dual_Conducting_Polymers

These contents are for calculating mobility of entities in an electric field

The main function is get_fit_window_plot.py
and the example file is: sample_MSD.txt

The file containing the MSD data should be formatted as in the sample file. MSDi is the MSD of entity i. t is the timestep of the given row. This can be in any desired unit.

To run the code type: "python get_fit_window_plot.py sample_MSD.txt nf E". Matplotlib, numpy, and scipy are all required, and the code is written using python 3 syntax.

Where ns is the conversion factor for the t column into femtosecons (1,000 for picosenconds, 1,000,000 for nansecond, etc.). E is the strength of the applied field in GigaVolts/m.

Currently, when you run the code, the data will be read and presented to you in batches of 5 MSD curves on one graph, along with the equation fit for each of those entities. Upon inspection, and which you decide can be fit well, close the window. The code will prompt you to type how many of those curves to fit; enter the value and hit enter. Next you will type which curves you would like to use. an entry begins with the curve number (taken from the legend in the plot), and if you wish to use the entire data of that curve hit enter and the entry for that curve is complete. An example screen would look like:

How many do we use?
 2
1
fit=  [890740898.29667747] 2.98452826808e-06
4
fit=  [39168056.593792528] 6.25843883039e-07

where we wanted to use curves 1 and 4 so we type we want to use 2 curves, then enter curves 1 and 4


If you only wish to fit a select portion of the curve you can add the randge of time (in nanoseconds) you want to fit. using the same example as above but fitting to a time window:

How many do we use?
 2
1 0 3
fit=  [890740898.29667747] 2.98452826808e-06
4 0 5
fit=  [39168056.593792528] 6.25843883039e-07

where we now we fit curves 1 and 4 to the time intervals (0,3) and (0,5), respectively. Once you type in as many entries as you said you wanted to use, the code will recalculate the fits, and replot the data with all the original MSD curves along with the new fit curves. once you have inspected the curves and approve of how they are fit, close the window and the code will ask if the fit is satisfactory. If yes, type 1 and hit enter. If no, type 0, hit enter and the code will replot all the original fits and you can adjust which curves you want to fit, and what time intervals to use.

If the fit is satisfactory the code will proceed to the next batch of curves, and the process repeats until all curves have been examined and satisfactorly fit. Once done, the code will take the fits from the selected curves, calculate the drift velocity, and calculate and average mobility in the field in units of m^2/Vs.
