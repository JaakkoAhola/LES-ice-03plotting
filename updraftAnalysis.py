#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 10:29:56 2020

@author: Jaakko Ahola, Finnish Meteorological Institute
@licence: MIT licence Copyright
"""


import matplotlib
import numpy
import pandas
import time
import sys
import os
import pathlib
from matplotlib.patches import Patch
import xarray

sys.path.append("../LES-03plotting")
from Figure import Figure
from Plot import Plot
from PlotTweak import PlotTweak
from Colorful import Colorful
from Data import Data
from SimulationDataAnalysis import SimulationDataAnalysis
from Simulation import Simulation
matplotlib.use("PDF")

class ManuscriptFigures:
    
    def __init__(self, figurefolder, datafolder, xend = 33.0):        
        
        self.figurefolder = figurefolder
        self.datafolder = pathlib.Path(datafolder)
        self.xstart = 2.1
        self.xend = xend
        
        self.simulation = Simulation( self.datafolder,"Prognostic", Colorful.getDistinctColorList("red"))

        print("Folders", self.figurefolder, self.datafolder)
        self.simulation.setAUXDataset("AeroB", "case_isdac_LVL5_3D_iceD_inter_48h_S_Nabb.nc")
        self.simulation.setAUXDataset("CloudB", "case_isdac_LVL5_3D_iceD_inter_48h_S_Ncbb.nc")
        self.simulation.setAUXDataset("IceB", "case_isdac_LVL5_3D_iceD_inter_48h_S_Nibb.nc")
        self.simulation.setAUXDataset("Updraft", "case_isdac_LVL5_3D_iceD_inter_48h_S_w.nc")


        self.simulation.setTimeCoordToHours()
        self.simulation.sliceByTimeAUXDataset(self.xstart, self.xend)






    def figureUpdraft(self):
        packing = 7


        yend = 1.5


        aeroAnalysis  = SimulationDataAnalysis( self.simulation, "S_Nabb", "AeroB")
        cloudAnalysis = SimulationDataAnalysis( self.simulation, "S_Ncbb", "CloudB")
        iceAnalysis   = SimulationDataAnalysis( self.simulation, "S_Nibb","IceB")


        aeroAnalysis.renameAUXCoordSizeBinB()
        cloudAnalysis.renameAUXCoordSizeBinB()
        iceAnalysis.renameAUXCoordSizeBinB()

        if packing <7:
            aeroAnalysis.packFilteredAUXVariablewithSizeBinCoords(packing)
            cloudAnalysis.packFilteredAUXVariablewithSizeBinCoords(packing)
            iceAnalysis.packFilteredAUXVariablewithSizeBinCoords(packing)

        aero = aeroAnalysis.simulation.AUXDatasets["AeroB"]["S_Nabb"]
        cloud = cloudAnalysis.simulation.AUXDatasets["CloudB"]["S_Ncbb"]
        ice = iceAnalysis.simulation.AUXDatasets["IceB"]["S_Nibb"]

        updraft = self.simulation.AUXDatasets["Updraft"]["w"]

        updraft = updraft.rename({"ym":"yt"})
        updraft.yt.values = updraft.yt.values-25.

        dataset = xarray.merge([aero,cloud,ice,updraft])
        #print(dataset)
        # print(" ")
        # print(aero)
        # print(" ")
        # print(cloud)
        # print(" ")
        # print(ice)
        # print(" ")

        # print(aero.sel(zt=200, method = "nearest"))



        aeroColor = Colorful.getDistinctColorList("red")
        cloudColor = Colorful.getDistinctColorList("navy")
        iceColor = Colorful.getDistinctColorList("cyan")

        yticks = [0, 0.5, 1, 1.5]


        heightList = [0,400,800]

        draftLimit = 1e-4

        for height in heightList:

            fig = Figure(self.figurefolder,"figureUpdraft" + str(height), ncols = 2, nrows =packing, figsize = [7,14])
            print("figsize", fig.getFigSize())
            datasetHeight = dataset.sel(zt = height, method="nearest")
            
            for draftIndex in range(1):
                if draftIndex == 0:
                    dataDraft = datasetHeight.where(datasetHeight["w"] > draftLimit, drop=True)
                    drafType = "Up-draft"
                else:
                    dataDraft = datasetHeight.where(datasetHeight["w"] < -draftLimit, drop=True)
                    drafType = "Down-draft"

                for bini in range(packing):
                    print("")
                    print("height", height, drafType,"bini", bini)

                    ax = fig.getAxes( bini*2 + draftIndex)

                    aeroHeight = dataDraft["S_Nabb"].mean(dim=["xt","yt"], skipna = True)
                    cloudHeight = dataDraft["S_Ncbb"].mean(dim=["xt","yt"], skipna = True)
                    iceHeight  = dataDraft["S_Nibb"].mean(dim=["xt","yt"], skipna = True)

                    aeroBin = aeroHeight[:,bini]
                    cloudBin = cloudHeight[:,bini]
                    iceBin = iceHeight[:,bini]

                    print(aeroBin)
                    print("")

                    print(cloudBin)
                    print("")

                    print(iceBin)
                    print("")

                    totalBin = aeroBin + cloudBin + iceBin

                    aeroFrac = aeroBin/totalBin
                    cloudFrac = cloudBin/totalBin
                    iceFrac = iceBin/totalBin

                    pointZero = totalBin.sel(time = self.xstart, method ="nearest")

                    pointEnd = totalBin.sel(time = self.xend, method ="nearest").values

                    totalBinRelative  = totalBin / pointZero


                    aeroFrac.plot(ax=ax, color = aeroColor)
                    cloudFrac.plot(ax=ax, color = cloudColor)
                    iceFrac.plot(ax=ax, color = iceColor)
                    totalBinRelative.plot(ax = ax, color = "black")

                    # if bini == (packing - 1):
                    #     bininame = str(bini + 1 ) + " - 7"
                    # else:
                    #     bininame = str(bini +1)

                    bininame = str(bini +1)
                    if True:
                        print("pointZero",pointZero, type(pointZero))
                        print("pointEnd", pointEnd, type(pointEnd))
                        label = " ".join([drafType, "Bin", bininame + ",", "Total", r"$N_0$",  str(int(pointZero)) + ",", "\nMin", r"$N$", str(int(pointEnd)), "$(kg^{-1})$"  ])
                        facecolor = "black"

                        PlotTweak.setArtist(ax, {label:facecolor}, loc = (0.01, 0.74), framealpha = 0.8)

                        if bini == 0:
                            collectionOfLabelsColors = {"Aerosol": aeroColor, "Cloud": cloudColor, "Ice": iceColor}
                            PlotTweak.setArtist(ax, collectionOfLabelsColors, ncol = 3, loc = (0.47,1.02))

                    ##############
                    ax.set_title("")
                    PlotTweak.setXLim(ax, start = self.xstart, end = self.xend)
                    PlotTweak.setYLim(ax, end = yend)

                    PlotTweak.setYticks(ax, yticks)
                    yShownLabelsBoolean = PlotTweak.setYLabels(ax, yticks, end = 1, interval = 1, integer=False)
                    PlotTweak.setYTickSizes(ax, yShownLabelsBoolean)

                    xticks = PlotTweak.setXticks(ax, start = self.xstart, end = self.xend, interval = 1)
                    shownLabelsBoolean = PlotTweak.setXLabels(ax, xticks, end = self.xend, interval = 4)
                    PlotTweak.setXTickSizes(ax, shownLabelsBoolean)

                    PlotTweak.setYaxisLabel(ax,"")
                    # if bini in [2,3]:
                    #     PlotTweak.setXaxisLabel(ax,"Time", "h")
                    # else:
                    #     PlotTweak.setXaxisLabel(ax,"")
                    #     PlotTweak.hideXTickLabels(ax)

                    # if bini in [1,3]:
                    #     PlotTweak.hideYTickLabels(ax)
                    ###########

                # end bini for loop
            # end draftIndex for loop
            fig.save()
        # end height for loop





def main(folder = os.environ["SIMULATIONFIGUREFOLDER"], datafolder = "/home/aholaj/Data/BinnedData", xend = 33.0):

    figObject = ManuscriptFigures(folder, datafolder, xend)


    if True:
        figObject.figureUpdraft()

if __name__ == "__main__":
    start = time.time()
    try:
        main( folder = sys.argv[1], datafolder = sys.argv[2], xend= sys.argv[3])
    except IndexError:
        main()
    
    end = time.time()
    print("Script completed in " + str(round((end - start),0)) + " seconds")
