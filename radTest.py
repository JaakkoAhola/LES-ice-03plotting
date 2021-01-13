#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 10:56:04 2020

@author: Jaakko Ahola, Finnish Meteorological Institute
@licence: MIT licence Copyright


"""
import matplotlib
import numpy
import pandas
import time
import sys
import os
from matplotlib.patches import Patch

sys.path.append(os.environ["LESMAINSCRIPTS"])
from InputSimulation import InputSimulation
from Figure import Figure
from Plot import Plot
from PlotTweak import PlotTweak
from Colorful import Colorful
from Data import Data
from SimulationDataAnalysis import SimulationDataAnalysis

from manuscriptFigures import ManuscriptFigures

class FigureRadSimulationList():

    def __init__(self):
        self.allSimulations = []
        self.ice0simulations = []
        self.ice1simulations = []
        self.ice4simulations = []
        models = ["COSMO","DHARMA-2M","DHARMA-bin","METO","RAMS", "SAM-2M", "SAM-bin", "UCLALES", "UCLALES-SB","WRFLES","WRFLES-PSU"]

        self.uclalesSimulations = ["ICE0_8h", "ICE1_8h", "ICE4_8h", "RadDamp"]


        for i in models:
            self.ice0simulations.append( i + "_" + "ice0")
            self.ice1simulations.append( i + "_" + "ice1")
            self.ice4simulations.append( i + "_" + "ice4")
        self.ice0simulations.append( self.uclalesSimulations[0] )
        self.ice1simulations.append( self.uclalesSimulations[1] )
        self.ice4simulations.append( self.uclalesSimulations[2] )
        self.ice4simulations.append( self.uclalesSimulations[3] )

        self.allSimulations = self.ice0simulations + self.ice1simulations + self.ice4simulations


    def getIce0Simulations(self):
        return self.ice0simulations
    def getIce1Simulations(self):
        return self.ice1simulations
    def getIce4Simulations(self):
        return self.ice4simulations

    def getAllSimulations(self):
        return self.allSimulations

    def getUCLALESSALSASimulations(self):
        return self.uclalesSimulations

figObject = ManuscriptFigures("/home/aholaj/Nextcloud/figures_updated/manuscriptSimulationData_Rad.csv",
                              os.environ["SIMULATIONFIGUREFOLDER"])

fig = Figure(figObject.figurefolder,"figure2RAD", ncols = 2, nrows = 3)

simulationList = FigureRadSimulationList()


for k in simulationList.getAllSimulations():
    try:
        figObject.simulationCollection[k].getTSDataset()
        figObject.simulationCollection[k].setTimeCoordToHours()
    except FileNotFoundError:
        if "ovchinnikov" in str(figObject.simulationCollection[k].getFolder()).lower():
            print("Ovchinnikov data is not available. Continue with existingsimulations")
            continue
        else:
            print("{0} data missing from {1}, make sure you have set your\
                      folder and environment variables correctly".format(k, figObject.simulationCollection[k].getFolder()))

    if figObject.simulationCollection[k].getLabel() == "BULK":
        figObject.simulationCollection[k].setZorder(1)
    elif figObject.simulationCollection[k].getLabel() == "BIN":
        figObject.simulationCollection[k].setZorder(2)
    else:
        figObject.simulationCollection[k].setZorder(3)

for k in simulationList.getUCLALESSALSASimulations():
    figObject.simulationCollection[k].setLineWidth(matplotlib.rcParams["lines.linewidth"]*1.5)
    figObject.simulationCollection[k].setColor(Colorful.getDistinctColorList("green"))

figObject.simulationCollection["RadDamp"].setColor(Colorful.getDistinctColorList("red"))
lwpYlimit = 60
iwpYlimit = 21
yPositionFrac = 0.91

# figA
ax = fig.getAxes(0)
for k in simulationList.getIce0Simulations():
    Plot.getTimeseries(ax,
                       figObject.simulationCollection[k],
                       "lwp_bar", conversionFactor=1000.)
PlotTweak.setAnnotation(ax, "a) ICE0 liquid water path", xPosition=0.2, yPosition= lwpYlimit*yPositionFrac)

# figB
ax = fig.getAxes(2)
for k in simulationList.getIce1Simulations():
    Plot.getTimeseries(ax,
                       figObject.simulationCollection[k],
                       "lwp_bar", conversionFactor=1000.)
PlotTweak.setAnnotation(ax, "b) ICE1 liquid water path", xPosition=0.2, yPosition= lwpYlimit*yPositionFrac)

# figC
ax = fig.getAxes(3)
for k in simulationList.getIce1Simulations():
    Plot.getTimeseries(ax,
                       figObject.simulationCollection[k],
                       "iwp_bar", conversionFactor=1000.)
PlotTweak.setAnnotation(ax, "c) ICE1 ice water path", xPosition=0.2, yPosition= iwpYlimit*yPositionFrac)

# figD
ax = fig.getAxes(4)
for k in simulationList.getIce4Simulations():
    Plot.getTimeseries(ax,
                       figObject.simulationCollection[k],
                       "lwp_bar", conversionFactor=1000.)

PlotTweak.setAnnotation(ax, "d) ICE4 liquid water path", xPosition=0.2, yPosition= lwpYlimit*yPositionFrac)

# figE
ax = fig.getAxes(5)
for k in simulationList.getIce4Simulations():
    Plot.getTimeseries(ax,
                       figObject.simulationCollection[k],
                       "iwp_bar", conversionFactor=1000.)

PlotTweak.setAnnotation(ax, "e) ICE4 ice water path", xPosition=0.2, yPosition= iwpYlimit*yPositionFrac)

for i in [0,2,4]: #LWP figures
    ax = fig.getAxes(i)
    end = lwpYlimit
    PlotTweak.setYaxisLabel(ax,"LWP", "g\ m^{-2}")
    PlotTweak.setYLim(ax, end = end)
    ticks = PlotTweak.setYticks(ax, end = end, interval = 2)
    shownLabelsBoolean = PlotTweak.setYLabels(ax, ticks, end = end, interval = 10)
    PlotTweak.setYTickSizes(ax, shownLabelsBoolean)

for i in [1,3,5]: #IWP figures
    ax = fig.getAxes(i)
    PlotTweak.setYaxisLabel(ax,"IWP", "g\ m^{-2}")
    end = iwpYlimit
    PlotTweak.setYLim(ax, end = end)
    ticks = PlotTweak.setYticks(ax, end = end, interval = 1)
    shownLabelsBoolean = PlotTweak.setYLabels(ax, ticks, end = end, interval = 3)
    PlotTweak.setYTickSizes(ax, shownLabelsBoolean)

for i in [0,2,3,4,5]: #all timeseries figures
    ax = fig.getAxes(i)
    end = 8
    PlotTweak.setXLim(ax, end = end)
    Plot.getVerticalLine(ax, 2)

    if i in [4,5]:
        PlotTweak.setXaxisLabel(ax,"Time", "h")
    else:
         PlotTweak.setXaxisLabel(ax,"")
         PlotTweak.hideXTickLabels(ax)

    ticks = PlotTweak.setXticks(ax, end = end, interval = 0.5, integer=False)
    ticks = [int(k) for k in  ticks]
    shownLabelsBoolean = PlotTweak.setXLabels(ax, ticks, end = end, interval = 2)
    PlotTweak.setXTickSizes(ax, shownLabelsBoolean)

# empty space
ax = fig.getAxes(1)
ax.axis("off")
legend_elements = [Patch(facecolor=figObject.simulationCollection["ICE0_8h"].getColor(),
                     label='UCLALES-SALSA'),
               Patch(facecolor=figObject.simulationCollection["SAM-bin_ice0"].getColor(),
                     label='BIN'),
               Patch(facecolor=figObject.simulationCollection["COSMO_ice0"].getColor(),
                     label='BULK'),
               Patch(facecolor=figObject.simulationCollection["RadDamp"].getColor(),
                     label='RadDamping')]

ax.legend(handles=legend_elements, loc='center', frameon = True, framealpha = 1.0)

fig.save()
