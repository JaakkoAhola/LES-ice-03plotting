#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 13:36:10 2020

@author: Jaakko Ahola, Finnish Meteorological Institute
@licence: MIT licence Copyright
"""

import matplotlib
import numpy
import pandas
import time
import sys
from matplotlib.patches import Patch
from matplotlib.ticker import ScalarFormatter,AutoMinorLocator

sys.path.append("/home/aholaj/Nextcloud/000_WORK/000_Codex/LES-superfolder/LES-03plotting")
from InputSimulation import InputSimulation
from Figure import Figure
from Plot import Plot
from PlotTweak import PlotTweak
from Colorful import Colorful
from Data import Data
from SimulationDataAnalysis import SimulationDataAnalysis

class ManuscriptFigures:
    
    def __init__(self, simulationDataFrameCSVFile, figurefolder):
        simulationDataFrame = pandas.read_csv(simulationDataFrameCSVFile)

        # set simulation data as dictionary
        self.simulationCollection = InputSimulation.getSimulationCollection( simulationDataFrame )
        self.figurefolder = figurefolder

    def figure2(self):
        
        # create figure object
        fig = Figure(self.figurefolder,"figure2", ncols = 2, nrows = 3)
        
        simulationList = Figure2SimulationList()
        
        
        
        
        for k in simulationList.getAllSimulations():
            self.simulationCollection[k].getTSDataset()
            self.simulationCollection[k].setTimeCoordToHours()
        
        for k in simulationList.getUCLALESSALSASimulations():
            self.simulationCollection[k].setLineWidth(matplotlib.rcParams["lines.linewidth"]*2)
            self.simulationCollection[k].setColor(Colorful.getDistinctColorList("green"))
        
        lwpYlimit = 60
        iwpYlimit = 21
        yPositionFrac = 0.92
        # figA
        ax = fig.getAxes(0)
        for k in simulationList.getIce0Simulations():
            Plot.getTimeseries(ax,
                               self.simulationCollection[k],
                               "lwp_bar", conversionFactor=1000.)
        PlotTweak.setAnnotation(ax, "a) ICE0 liquid water path", xPosition=0.2, yPosition= lwpYlimit*yPositionFrac)
        
        # figB
        ax = fig.getAxes(2)
        for k in simulationList.getIce1Simulations():
            Plot.getTimeseries(ax,
                               self.simulationCollection[k],
                               "lwp_bar", conversionFactor=1000.)
        PlotTweak.setAnnotation(ax, "b) ICE1 liquid water path", xPosition=0.2, yPosition= lwpYlimit*yPositionFrac)
            
        
        # figC
        ax = fig.getAxes(3)
        for k in simulationList.getIce1Simulations():
            Plot.getTimeseries(ax,
                               self.simulationCollection[k],
                               "iwp_bar", conversionFactor=1000.)
        PlotTweak.setAnnotation(ax, "c) ICE1 ice water path", xPosition=0.2, yPosition= iwpYlimit*yPositionFrac)
            
            
        
        # figD
        ax = fig.getAxes(4)
        for k in simulationList.getIce4Simulations():
            Plot.getTimeseries(ax,
                               self.simulationCollection[k],
                               "lwp_bar", conversionFactor=1000.)
        PlotTweak.setAnnotation(ax, "d) ICE1 liquid water path", xPosition=0.2, yPosition= lwpYlimit*yPositionFrac)
        
        # figE
        ax = fig.getAxes(5)
        for k in simulationList.getIce4Simulations():
            Plot.getTimeseries(ax,
                               self.simulationCollection[k],
                               "iwp_bar", conversionFactor=1000.)
        PlotTweak.setAnnotation(ax, "e) ICE1 ice water path", xPosition=0.2, yPosition= iwpYlimit*yPositionFrac)
        
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
        
            # set xticks
            ticks = PlotTweak.setXticks(ax, end = end, interval = 0.5, integer=False)
            # set xlabels
            ticks = [int(k) for k in  ticks]
            shownLabelsBoolean = PlotTweak.setXLabels(ax, ticks, end = end, interval = 2)
            # set xtick sizes
            PlotTweak.setXTickSizes(ax, shownLabelsBoolean)
            
        
        # empty space
        ax = fig.getAxes(1)
        ax.axis("off")
        legend_elements = [Patch(facecolor=self.simulationCollection["ICE0_8h"].getColor(),
                             label='UCLALES-SALSA'),
                       Patch(facecolor=self.simulationCollection["COSMO_ice0"].getColor(),
                             label='BIN'),
                       Patch(facecolor=self.simulationCollection["SAM-bin_ice0"].getColor(),
                             label='BULK')]
    
        ax.legend(handles=legend_elements, loc='center', frameon = True, framealpha = 1.0)
            
            
        fig.save()
        
    def figure3(self):
        
        # create figure object
        fig = Figure(self.figurefolder,"figure3", ncols = 2, nrows = 1, figsize = [12/2.54, 7/2.54],bottom = 0.18, top=0.8)
        
        simulationList = []
        for i in range(1,7):
            case = str(i)
            pretext = "ICE" + case 
            if case in ["5", "6"]:
                posttext = "8h"
            else:
                posttext = "24h"
            simulationList.append(pretext+ "_" + posttext)
        
        for k in simulationList:
            self.simulationCollection[k].getTSDataset()
            self.simulationCollection[k].setTimeCoordToHours()
        
        
        
        lwpYlimit = 60
        iwpYlimit = 22
        
        xPosition = 0.5
        yPositionFrac = 0.92
        
        # figA
        
        ax = fig.getAxes(0)
        for k in simulationList:
            Plot.getTimeseries(ax,
                           self.simulationCollection[k],
                           "lwp_bar", conversionFactor=1000.)
        end = lwpYlimit
        PlotTweak.setAnnotation(ax, "a) Liquid water path", xPosition=xPosition, yPosition= yPositionFrac*end)
        PlotTweak.setYaxisLabel(ax,"LWP", "g\ m^{-2}")
        PlotTweak.setYLim(ax, end = end)
        PlotTweak.setYLim(ax, end = end)  
        ticks = PlotTweak.setYticks(ax, end = end, interval = 2)
        shownLabelsBoolean = PlotTweak.setYLabels(ax, ticks, end = end, interval = 10)
        PlotTweak.setYTickSizes(ax, shownLabelsBoolean)
        
        
        # figB
        ax = fig.getAxes(1)
        for k in simulationList:
            Plot.getTimeseries(ax,
                               self.simulationCollection[k],
                               "iwp_bar", conversionFactor=1000.)
        end = iwpYlimit
        PlotTweak.setAnnotation(ax, "b) Ice water path", xPosition=5, yPosition= yPositionFrac*end)
        PlotTweak.setYaxisLabel(ax,"IWP", "g\ m^{-2}")
        PlotTweak.setYLim(ax, end = end)
        end = iwpYlimit
        PlotTweak.setYLim(ax, end = end)  
        ticks = PlotTweak.setYticks(ax, end = end, interval = 1)
        shownLabelsBoolean = PlotTweak.setYLabels(ax, ticks, end = end, interval = 3)
        PlotTweak.setYTickSizes(ax, shownLabelsBoolean)
    
        for i in [0,1]:
            ax = fig.getAxes(i)
            end = 24
            PlotTweak.setXLim(ax, end = end)
            Plot.getVerticalLine(ax, 2)
            Plot.getVerticalLine(ax, 8)
            PlotTweak.setXaxisLabel(ax,"Time", "h")
            
            # set xticks
            ticks = PlotTweak.setXticks(ax, end = end, interval = 1, integer=False)
            # set xlabels
            ticks = [int(k) for k in  ticks]
            shownLabelsBoolean = PlotTweak.setXLabels(ax, ticks, end = end, interval = 4)
            # set xtick sizes
            PlotTweak.setXTickSizes(ax, shownLabelsBoolean)                
        
        ax = fig.getAxes(0)
        legend_elements = []
        
        for k in simulationList:
            legend_elements.append(Patch(facecolor=self.simulationCollection[k].getColor(),
                             label=self.simulationCollection[k].getLabel()))
                       
    
        ax.legend(handles=legend_elements, loc=(0.5,1.05), frameon = True, framealpha = 1.0, ncol=3)
            
        fig.save()
        
        
        
    def figure4(self):
        
        fig = Figure(self.figurefolder,"figure4", ncols = 2, nrows = 3, left=0.15, wspace=0.4)
        
        simulationList = ["Prognostic_48h", "ICE4_24h"]
        
        
        for k in simulationList:
            self.simulationCollection[k].getTSDataset()
            self.simulationCollection[k].getNCDataset()
            self.simulationCollection[k].setTimeCoordToHours()
        
        
        
        lwpYlimit = 60
        iwpYlimit = 21
        
        xPosition = 0.52
        yPositionFrac = 0.91
        
        # figA
        
        ax = fig.getAxes(0)
        for k in simulationList:
            Plot.getTimeseries(ax,
                           self.simulationCollection[k],
                           "lwp_bar", conversionFactor=1000.)
        end = lwpYlimit
        PlotTweak.setAnnotation(ax, "a) Liquid water path", xPosition=xPosition, yPosition= yPositionFrac*end)
        PlotTweak.setYaxisLabel(ax,"LWP", "g\ m^{-2}")
        PlotTweak.setYLim(ax, end = end)
        PlotTweak.setYLim(ax, end = end)  
        ticks = PlotTweak.setYticks(ax, end = end, interval = 2)
        shownLabelsBoolean = PlotTweak.setYLabels(ax, ticks, end = end, interval = 10)
        PlotTweak.setYTickSizes(ax, shownLabelsBoolean)
        
        
        # figB
        ax = fig.getAxes(1)
        for k in simulationList:
            Plot.getTimeseries(ax,
                               self.simulationCollection[k],
                               "iwp_bar", conversionFactor=1000.)
        end = iwpYlimit
        PlotTweak.setAnnotation(ax, "b) Ice water path", xPosition=xPosition, yPosition= yPositionFrac*end)
        PlotTweak.setYaxisLabel(ax,"IWP", "g\ m^{-2}")
        PlotTweak.setYLim(ax, end = end)
        end = iwpYlimit
        PlotTweak.setYLim(ax, end = end)  
        ticks = PlotTweak.setYticks(ax, end = end, interval = 1)
        shownLabelsBoolean = PlotTweak.setYLabels(ax, ticks, end = end, interval = 3)
        PlotTweak.setYTickSizes(ax, shownLabelsBoolean)
            
            
        
        # figC
        ax = fig.getAxes(2)
        for k in simulationList:
            self.simulationCollection[k].getNCDataset()["Ni_ii_vol"].plot(ax = ax,
                                                                     color = self.simulationCollection[k].getColor(),
                                                                     label =  self.simulationCollection[k].getLabel(),
                                                                     linewidth = self.simulationCollection[k].getLineWidth())
        end = 5
        PlotTweak.setAnnotation(ax, "c) Ice number \nconcentration", xPosition=xPosition, yPosition= 4)
        PlotTweak.setYaxisLabel(ax,"N_{i}", "L^{-1}")
        PlotTweak.setYLim(ax, end = end)  
        ticks = PlotTweak.setYticks(ax, end = end, interval = 0.2, integer=False)
        
        shownLabelsBoolean = PlotTweak.setYLabels(ax, ticks, end = end, interval = 1)
        PlotTweak.setYTickSizes(ax, shownLabelsBoolean)
            
            
        # figD
        ax = fig.getAxes(3)
        for k in simulationList:
            Plot.getTimeseries(ax,
                               self.simulationCollection[k],
                               "Nc_ic", conversionFactor=1e-6)
        end = 170
        PlotTweak.setAnnotation(ax, "d) In-cloud CDNC", xPosition=xPosition, yPosition= yPositionFrac*end)
        PlotTweak.setYaxisLabel(ax,"CDNC", "mg^{-1}")
        PlotTweak.setYLim(ax, end = end)  
        ticks = PlotTweak.setYticks(ax, end = end, interval = 10)
        shownLabelsBoolean = PlotTweak.setYLabels(ax, ticks, end = end, interval = 50)
        PlotTweak.setYTickSizes(ax, shownLabelsBoolean)
            
        # figE
        ax = fig.getAxes(4)
        for k in simulationList:
            for muuttuja in ["zb", "zc"]:
                self.simulationCollection[k].getTSDataset()[muuttuja].where(self.simulationCollection[k].getTSDataset()[muuttuja] > 0).plot(ax = ax,
                                                                     color = self.simulationCollection[k].getColor(),
                                                                     label =  self.simulationCollection[k].getLabel(),
                                                                     linewidth = self.simulationCollection[k].getLineWidth())
        end = 1000
        PlotTweak.setAnnotation(ax, "e) Cloud top and base", xPosition=xPosition, yPosition= yPositionFrac*end)
        PlotTweak.setYaxisLabel(ax,"Height", "m")
        PlotTweak.setYLim(ax, end = end)    
        ticks = PlotTweak.setYticks(ax, end = end, interval = 100)
        shownLabelsBoolean = PlotTweak.setYLabels(ax, ticks, end = end, interval = 200)
        PlotTweak.setYTickSizes(ax, shownLabelsBoolean)
                    
        for i in range(5):
            ax = fig.getAxes(i)
            end = 32
            PlotTweak.setXLim(ax, end = end)
            Plot.getVerticalLine(ax, 2)
            Plot.getVerticalLine(ax, 24)
            
            
            if i in [3,4]:
                PlotTweak.setXaxisLabel(ax,"Time", "h")
            else:
                PlotTweak.setXaxisLabel(ax,"")
                PlotTweak.hideXTickLabels(ax)
        
            ticks = PlotTweak.setXticks(ax, end = end, interval = 0.5)
            shownLabelsBoolean = PlotTweak.setXLabels(ax, ticks, end = end, interval = 4)
            PlotTweak.setXTickSizes(ax, shownLabelsBoolean)
            
        
        # empty space
        ax = fig.getAxes(5)
        ax.axis("off")
        legend_elements = [Patch(facecolor=self.simulationCollection["Prognostic_48h"].getColor(),
                             label='Prognostic ice'),
                       Patch(facecolor=self.simulationCollection["ICE4_24h"].getColor(),
                             label='ICE4')]
    
        ax.legend(handles=legend_elements, loc='center', frameon = True, framealpha = 1.0)
            
        fig.save()
        
    
    def figure5(self):
        
        fig = Figure(self.figurefolder,"figure5", ncols = 2, nrows = 2, wspace=0.1, bottom = 0.14, left=0.15, top=0.98)
        
        simulation = "Prognostic_48h"
        
        self.simulationCollection[simulation].getPSDataset()
        self.simulationCollection[simulation].getTSDataset()
        self.simulationCollection[simulation].setTimeCoordToHours()
            
        logaritmicLevels = PlotTweak.getLogaritmicTicks(-17,-9, includeFives = True)
        end = 1000
        yPositionFrac = 0.9
        xPosition = 0.5
        orange = Colorful.getDistinctColorList("orange")
        
        annotation = ["a) Dust in aerosols", "b) Dust in cloud droplets", "c) Dust in ice crystals"]
        for ind, muuttuja in enumerate(["P_cDUa", "P_cDUc", "P_cDUi"]):
            ax = fig.getAxes(ind)
            data = self.simulationCollection[simulation].getPSDataset()[muuttuja]
            data.values  = numpy.log10(data.values)
            
            im = data.plot.contourf("time","zt", ax = ax, levels=logaritmicLevels, add_colorbar = False)
            self.simulationCollection[simulation].getTSDataset()["zb"].plot(ax = ax,
                                                                     color = orange,
                                                                     linewidth = self.simulationCollection[simulation].getLineWidth())
            
            ax = Plot.getContourLine(ax, self.simulationCollection[simulation], "P_RHi", 100)
            PlotTweak.setAnnotation(ax, annotation[ind], xPosition=xPosition, yPosition= yPositionFrac*end)
            
            end = 1000
            PlotTweak.setYLim(ax, end = end)    
            ticks = PlotTweak.setYticks(ax, end = end, interval = 100)
            shownLabelsBoolean = PlotTweak.setYLabels(ax, ticks, end = end, interval = 200)
            PlotTweak.setYTickSizes(ax, shownLabelsBoolean)
            
            
            
            
        
        
        for i in range(3):
            ax = fig.getAxes(i)
            end = 33.05
            PlotTweak.setXLim(ax, end = end)
            Plot.getVerticalLine(ax, 2)
            
            if i in [1,2]:
                PlotTweak.setXaxisLabel(ax,"Time", "h")
            else:
                PlotTweak.setXaxisLabel(ax,"")
                PlotTweak.hideXTickLabels(ax)
            
            if i == 1:
                PlotTweak.hideYTickLabels(ax)
            
            if i in [0,2]:
                PlotTweak.setYaxisLabel(ax,"Height", "m")
            else:
                PlotTweak.setYaxisLabel(ax,"")
                PlotTweak.hideYTickLabels(ax)
        
            ticks = PlotTweak.setXticks(ax, end = end, interval = 1)
            shownLabelsBoolean = PlotTweak.setXLabels(ax, ticks, end = end, interval = 4)
            PlotTweak.setXTickSizes(ax, shownLabelsBoolean)
            
        
        # empty space
        ax = fig.getAxes(3)
        ax.axis("off")
        legend_elements = [Patch(facecolor="black",
                             label='RH over ice 100%'),
                       Patch(facecolor=orange,
                             label='Cloud base')]
    
        ax.legend(handles=legend_elements, loc=(0.07,0), frameon = True, framealpha = 1.0)
        cax = matplotlib.pyplot.axes([0.6, 0.4, 0.33, 0.03])
        Plot.getColorBar(im,cax,logaritmicLevels)
        
        shownLabelsBoolean = Data.getMaskedList(logaritmicLevels, numpy.arange(-11,-19,-1))
        PlotTweak._hideLabels(cax.xaxis, shownLabelsBoolean)
        PlotTweak.setXTickSizes(cax, shownLabelsBoolean)
        fig.save()
    
    def figure6(self):
        packing = 4
        
        simulation = "Prognostic_48h"
        simulCol = self.simulationCollection[simulation]
        simulCol.getPSDataset()
        
        
        simulCol.setTimeCoordToHours()
        
        fig = Figure(self.figurefolder,"figure6", ncols = 2, nrows = 2, hspace=0.1, bottom = 0.14, left=0.05, top=0.9, wspace = 0.06)
        
        aeroAnalysis  = SimulationDataAnalysis( simulCol, "P_Nabb")
        cloudAnalysis = SimulationDataAnalysis( simulCol, "P_Ncbb")
        iceAnalysis   = SimulationDataAnalysis( simulCol, "P_Nibb")
    
        
        aeroAnalysis.filterPSVariableInCloud()
        cloudAnalysis.filterPSVariableInCloud()
        iceAnalysis.filterPSVariableInCloud()
        
        aeroAnalysis.renamePSCoordSizeBinB()
        cloudAnalysis.renamePSCoordSizeBinB()
        iceAnalysis.renamePSCoordSizeBinB()
        
        
        aeroAnalysis.packFilteredPSVariablewithSizeBinCoords(packing)
        cloudAnalysis.packFilteredPSVariablewithSizeBinCoords(packing)
        iceAnalysis.packFilteredPSVariablewithSizeBinCoords(packing)
        
        aero = simulCol.getPSDataset()[aeroAnalysis.getFilteredPackedVariableName()]
        cloud = simulCol.getPSDataset()[cloudAnalysis.getFilteredPackedVariableName()]
        ice = simulCol.getPSDataset()[iceAnalysis.getFilteredPackedVariableName()]
        
        total = aero + cloud + ice
        
        aeroColor = Colorful.getDistinctColorList("red")
        cloudColor = Colorful.getDistinctColorList("navy")
        iceColor = Colorful.getDistinctColorList("cyan")
        
        xstart = 2
        xend = 33.05
        
        yend = 1.5
        yticks = [0, 0.5, 1]
        
        print(total.values)
        
        figName = ["a)", "b)", "c)", "d)"]
        
        print("shape total", numpy.shape(total.values))
        for bini in range(packing):
            ax = fig.getAxes(bini)
            
            aeroBin = aero[:,bini]
            cloudBin = cloud[:,bini]
            iceBin = ice[:,bini]
            
            totalBin = total[:,bini]
            
            aeroFrac = aeroBin/totalBin
            cloudFrac = cloudBin/totalBin
            iceFrac = iceBin/totalBin
            print(" ")
            print(totalBin.sel(time = xstart, method ="nearest"))
            totalBinRelative  = totalBin / totalBin.sel(time = xstart, method ="nearest")
            
            
            
            aeroFrac.plot(ax=ax, color = aeroColor)
            cloudFrac.plot(ax=ax, color = cloudColor)
            iceFrac.plot(ax=ax, color = iceColor)
            print("totalBinRelative", max(totalBinRelative), min(totalBinRelative))
            totalBinRelative.plot(ax = ax, color = "black")
            
            if bini == (packing - 1):
                bininame = str(bini + 1 ) + " - 7"
            else:
                bininame = str(bini +1)
            
            if True:
                
                label = " ".join([figName[bini], "Bin", bininame + ",", "Total", r"$N_0$",  str(int(totalBin.values[0])) + ",", "\nMin", r"$N$", str(int(numpy.min(totalBin))), "$(kg^{-1})$"  ])
                facecolor = "black"
                
                PlotTweak.setArtist(ax, {label:facecolor}, loc = (0.01, 0.67), framealpha = 0.8)
            
                if bini == 0:
                    collectionOfLabelsColors = {"Aerosol": aeroColor, "Cloud": cloudColor, "Ice": iceColor}
                    PlotTweak.setArtist(ax, collectionOfLabelsColors, ncol = 3, loc = (0.5,1.05))
            
            ##############
            ax.set_title("")
            PlotTweak.setXLim(ax, start = xstart, end = xend)
            PlotTweak.setYLim(ax, end = yend)
            
            PlotTweak.setYticks(ax, yticks)
            yShownLabelsBoolean = PlotTweak.setYLabels(ax, yticks, end = 1, interval = 1, integer=False)
            PlotTweak.setYTickSizes(ax, yShownLabelsBoolean)
            
            xticks = PlotTweak.setXticks(ax, start = xstart, end = xend, interval = 1)
            shownLabelsBoolean = PlotTweak.setXLabels(ax, xticks, end = xend, interval = 4)
            PlotTweak.setXTickSizes(ax, shownLabelsBoolean)
            
            PlotTweak.setYaxisLabel(ax,"")
            if bini in [2,3]:
                PlotTweak.setXaxisLabel(ax,"Time", "h")
            else:
                PlotTweak.setXaxisLabel(ax,"")
                PlotTweak.hideXTickLabels(ax)
                
            if bini in [1,3]:
                PlotTweak.hideYTickLabels(ax)
            ###########
            
        fig.save()    
    
    def figure7(self):
        
        fig = Figure(self.figurefolder,"figure7", ncols = 2, nrows = 2, hspace=0.43, bottom = 0.14, left=0.15, top=0.98, wspace = 0.1)
        
        simulation = "Prognostic_48h"
        
        self.simulationCollection[simulation].getPSDataset()
        self.simulationCollection[simulation].setTimeCoordToHours()
            
        end = 1000
        yPositionFrac = 0.9
        xPositionFrac = 0.02
        
        annotation = ["a) Liquid water mixing ratio", "b) Ice mixing ratio", "c) Freezing rate"]
        muuttujalista = ["P_rl", "P_ri", "nucl_ni"]
        aika = [2,4,16,32]
        colorlist = ['#a6cee3','#1f78b4','#b2df8a','#33a02c'] #['#a63603', '#e6550d', '#fd8d3c', '#fdae6b'] #['#1b9e77','#d95f02','#7570b3','#e7298a']
        
        startList = [0, 0, -5]
        endList = [0.2, 0.025,2.5]
        
        
        xPositionList = []
        
        for ind, muuttuja in enumerate(muuttujalista):
            ax = fig.getAxes(ind)
            data = self.simulationCollection[simulation].getPSDataset()[muuttuja]
            
            xPositionList.append( numpy.abs(startList[ind]-endList[ind])*xPositionFrac + startList[ind])
            
            if muuttuja == "nucl_ni":
                data.values  = numpy.log10(data.values)
            else:
                data.values = data.values*1e3
                
            
            for aikaInd, aikapiste in enumerate(aika):
                
                data.sel(time = aikapiste, method="nearest").plot( y="zt", ax = ax, color = colorlist[aikaInd])
                
            
                
                
                
            #
    
        end = 1000
        for ind in range(3):
            ax = fig.getAxes(ind)
            ax.set_title("")
            PlotTweak.setYLim(ax, end = end)    
            ticks = PlotTweak.setYticks(ax, end = end, interval = 100)
            shownLabelsBoolean = PlotTweak.setYLabels(ax, ticks, end = end, interval = 200)
            PlotTweak.setYTickSizes(ax, shownLabelsBoolean)
            PlotTweak.setAnnotation(ax, annotation[ind], xPosition= xPositionList[ind], yPosition= yPositionFrac*end)
            
            
            
        
        # for i in range(3):
        #     if i in [0,2]:
        #         ax = fig.getAxes(i)
        #         PlotTweak.setYaxisLabel(ax,"Height", "m")
        #     else:
        #         PlotTweak.setYaxisLabel(ax,"")
        #         PlotTweak.hideYTickLabels(ax)
        
        
        ind = 0
        ax = fig.getAxes(ind)
        start = startList[ind]
        end = endList[ind]
        PlotTweak.setXLim(ax, start = ind, end = end) 
        PlotTweak.setXaxisLabel(ax,"", "g\ kg^{-1}")
        PlotTweak.setYaxisLabel(ax,"Height", "m")
        
        matplotlib.pyplot.setp(ax.get_xticklabels()[-1], visible=False)
        
        # ticks = PlotTweak.setXticks(ax, end = end, interval = 0.05/2, integer=False)
        # shownLabelsBoolean = PlotTweak.setXLabels(ax, ticks, end = end, interval = 0.05, integer = False)
        # PlotTweak.setXTickSizes(ax, shownLabelsBoolean)
        
        ind = 1
        ax = fig.getAxes(ind)
        start = startList[ind]
        end = endList[ind]
        PlotTweak.setXLim(ax, start = start, end = end) 
        PlotTweak.setXaxisLabel(ax,"", "g\ kg^{-1}")
        PlotTweak.setYaxisLabel(ax,"")
        PlotTweak.hideYTickLabels(ax)
        
        
        ticks = PlotTweak.setXticks(ax, end = end, interval = 0.005/2, integer=False)
        shownLabelsBoolean = PlotTweak.setXLabels(ax, ticks, end = end, interval = 0.005, integer = False)
        PlotTweak.setXTickSizes(ax, shownLabelsBoolean)
        matplotlib.pyplot.setp(ax.get_xticklabels()[-1], visible=False)
        
        
        ind = 2
        ax = fig.getAxes(ind)
        start = startList[ind]
        end = endList[ind]
        PlotTweak.setYaxisLabel(ax,"Height", "m")
        PlotTweak.setXaxisLabel(ax,"", "m^{-3}\ s^{-1}")
        PlotTweak.setXLim(ax, start = start, end = end) 
        
        xticksLog = numpy.arange(start, end + 0.1, 2.5)
        xlabelsLog = [r"$10^{" + "{0:.1f}".format(elem) + "}$" for elem in xticksLog]
        ax.set_xticks(xticksLog)
        ax.set_xticklabels(xlabelsLog)
    
        # ticks = PlotTweak.setXticks(ax, end = end, interval = 0.5)
        # shownLabelsBoolean = PlotTweak.setXLabels(ax, ticks, end = end, interval = 4)
        # PlotTweak.setXTickSizes(ax, shownLabelsBoolean)
        
        # empty space
        ax = fig.getAxes(3)
        ax.axis("off")
        legend_elements = [Patch(facecolor=colorlist[0], label='2 (h)'),
                       Patch(facecolor=colorlist[1], label='4 (h)'),
                       Patch(facecolor=colorlist[2], label='16 (h)'),
                       Patch(facecolor=colorlist[3], label='32 (h)')]
    
        ax.legend(handles=legend_elements, loc="center", frameon = True, framealpha = 1.0)
        fig.save()

class Figure2SimulationList():
    
    def __init__(self):
        self.allSimulations = []
        self.ice0simulations = []
        self.ice1simulations = []
        self.ice4simulations = []
        models = ["COSMO","DHARMA-2M","DHARMA-bin","METO","RAMS", "SAM-2M", "SAM-bin", "UCLALES", "UCLALES-SB","WRFLES","WRFLES-PSU"]
        
        self.uclalesSimulations = ["ICE0_8h", "ICE1_8h", "ICE4_8h"]
        
        
        for i in models:
            self.ice0simulations.append( i + "_" + "ice0")
            self.ice1simulations.append( i + "_" + "ice1")
            self.ice4simulations.append( i + "_" + "ice4")
        self.ice0simulations.append( self.uclalesSimulations[0] )
        self.ice1simulations.append( self.uclalesSimulations[1] )
        self.ice4simulations.append( self.uclalesSimulations[2] )
        
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
        
        
    
    

        

def main():

    figObject = ManuscriptFigures("/home/aholaj/Nextcloud/figures_updated/manuscriptSimulationData.csv", 
                                  "/home/aholaj/Nextcloud/figures_updated")
    
    if False:
        figObject.figure2()
    if False:
        figObject.figure3()
    if False:
        figObject.figure4()
    if False:
        figObject.figure5()
    if True:
        figObject.figure6()
    if False:
        figObject.figure7()        
    
    
     
    
if __name__ == "__main__":
    start = time.time()
    
    main()
    
    end = time.time()
    print("Script completed in " + str(round((end - start),0)) + " seconds")
