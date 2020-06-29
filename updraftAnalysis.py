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

    def __init__(self, figurefolder, datafolder, xstart = 2.1, xend = 33.0, draftLimit = 1e-3):

        self.figurefolder = pathlib.Path(figurefolder)
        self.figurefolder.mkdir( parents=True, exist_ok = True )

        self.datafolder = pathlib.Path(datafolder)
        self.xstart = xstart
        self.xend = xend

        self.draftLimit = draftLimit

        self.simulation = Simulation( self.datafolder,"Prognostic", Colorful.getDistinctColorList("red"))

        print("Folders", self.figurefolder, self.datafolder)
        self.simulation.setAUXDataset("AeroB", "case_isdac_LVL5_3D_iceD_inter_48h_S_Nabb.nc")
        self.simulation.setAUXDataset("CloudB", "case_isdac_LVL5_3D_iceD_inter_48h_S_Ncbb.nc")
        self.simulation.setAUXDataset("IceB", "case_isdac_LVL5_3D_iceD_inter_48h_S_Nibb.nc")
        self.simulation.setAUXDataset("Updraft", "case_isdac_LVL5_3D_iceD_inter_48h_S_w.nc")


        self.simulation.setTimeCoordToHours()
        self.simulation.sliceByTimeAUXDataset(self.xstart, self.xend)






    def figureUpdraftTimeseries(self):
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

        yticks = numpy.arange(0, 3.5+0.1, 0.5)


        heightList = {0:"Surface",
                200:"",
                400:"",
                355:"Always below cloud base",
                705:"Always in-cloud",
                785:"At the beginning in-cloud, At end above cloud top",
                850:"Always above cloud top"}

        for height in heightList:
            realHeight = dataset.zt.sel(zt = height, method="nearest").item()
            fig = Figure(self.figurefolder,"figureUpdraft_z_" + "{0:.0f}".format(realHeight),
                         ncols = 2, nrows =packing,
                         figsize = [8,16], wspace=0.06, left=0.05, bottom = 0.03, top=0.96, right=0.98)

            print("figsize", fig.getFigSize())
            datasetHeight = dataset.sel(zt = height, method="nearest")

            for draftIndex in range(2):
                if draftIndex == 0:
                    dataDraft = datasetHeight.where(datasetHeight["w"] > self.draftLimit, drop=True)
                    drafType = "Up-draft"
                else:
                    dataDraft = datasetHeight.where(datasetHeight["w"] < - self.draftLimit, drop=True)
                    drafType = "Down-draft"

                for bini in range(packing):
                    print("")
                    print("height", realHeight, drafType,"bini", bini)
                    axIndex =  bini*2 + draftIndex
                    ax = fig.getAxes(axIndex)

                    if dataDraft["w"].size == 0:
                        ax.axis("off")
                        continue

                    aeroHeight = dataDraft["S_Nabb"].mean(dim=["xt","yt"], skipna = True)
                    cloudHeight = dataDraft["S_Ncbb"].mean(dim=["xt","yt"], skipna = True)
                    iceHeight  = dataDraft["S_Nibb"].mean(dim=["xt","yt"], skipna = True)

                    aeroBin = aeroHeight[:,bini]
                    cloudBin = cloudHeight[:,bini]
                    iceBin = iceHeight[:,bini]

                    #print(aeroBin)
                    #print("")

                    #print(cloudBin)
                    #print("")

                    #print(iceBin)
                    #print("")

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
                        #print("pointZero",pointZero, type(pointZero))
                        #print("pointEnd", pointEnd, type(pointEnd))
                        label = " ".join([drafType, "Bin", bininame + ",", "Total", r"$N_0$",  str(int(pointZero)) + ",", "\nMin", r"$N$", str(int(pointEnd)), "$(kg^{-1})$"  ])
                        facecolor = "black"

                        PlotTweak.setArtist(ax, {label:facecolor}, loc = (0.01, 0.74), framealpha = 0.8)

                        if axIndex == 0:
                            collectionOfLabelsColors = {"Aerosol": aeroColor, "Cloud": cloudColor, "Ice": iceColor}
                            PlotTweak.setArtist(ax, collectionOfLabelsColors, ncol = 3, loc = (0.75,1.12))

                    ##############
                    ax.set_title("")
                    PlotTweak.setXLim(ax, start = self.xstart, end = self.xend)
                    PlotTweak.setYLim(ax, end = yend)

                    PlotTweak.setYticks(ax, yticks)
                    yShownLabelsBoolean = PlotTweak.setYLabels(ax, yticks, end = 3, interval = 1, integer=True)
                    PlotTweak.setYTickSizes(ax, yShownLabelsBoolean)

                    xticks = PlotTweak.setXticks(ax, start = self.xstart, end = self.xend, interval = 1)
                    shownLabelsBoolean = PlotTweak.setXLabels(ax, xticks, end = self.xend, interval = 4)
                    PlotTweak.setXTickSizes(ax, shownLabelsBoolean)

                    PlotTweak.setYaxisLabel(ax,"")
                    if axIndex in [12,13]:
                        PlotTweak.setXaxisLabel(ax,"Time", "h")
                    else:
                        PlotTweak.setXaxisLabel(ax,"")
                        PlotTweak.hideXTickLabels(ax)

                    if draftIndex == 1:
                        PlotTweak.hideYTickLabels(ax)
                    if axIndex == 0:
                        ax.text( 0.5*self.xend, yticks[-1]+yticks[1]*0.25, PlotTweak.getUnitLabel("Height\ " + f"{realHeight:.0f}", "m")+ " " +  heightList[height] + " limit: " + f"{self.draftLimit:.0e}"  , size=8)

                # end bini for loop
            # end draftIndex for loop
            fig.save()
        # end height for loop

    def figureUpdraftProfile(self):
        packing = 7


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

        aeroColor = Colorful.getDistinctColorList("red")
        cloudColor = Colorful.getDistinctColorList("navy")
        iceColor = Colorful.getDistinctColorList("cyan")

        xaxisend = 1.0
        yend = 1000.

        fig = Figure(self.figurefolder,"figureProfileDrafts",
                     ncols = 2, nrows =packing,
                     wspace=0.06, left=0.12, bottom = 0.03, top=0.96, right=0.98)
        print("figsize", fig.getFigSize())

        timeBegin = dataset.time.sel(time=self.xstart, method = "nearest").item()
        timeEnd = dataset.time.sel(time=self.xend, method = "nearest").item()
        datasetTime = dataset.sel(time = slice(timeBegin,timeEnd))

        for draftIndex in range(2):
            if draftIndex == 0:
                dataDraft = datasetTime.where(datasetTime["w"] > self.draftLimit, drop=True)
                drafType = "Up-draft"
            else:
                dataDraft = datasetTime.where(datasetTime["w"] < - self.draftLimit, drop=True)
                drafType = "Down-draft"

            for bini in range(packing):
                print("")
                print(drafType, "timeBegin", timeBegin, "timeEnd", timeEnd)
                axIndex =  bini*2 + draftIndex
                ax = fig.getAxes(axIndex)

                if dataDraft["w"].size == 0:
                    ax.axis("off")
                    continue

                aeroHeight = dataDraft["S_Nabb"].mean(dim=["xt","yt", "time"], skipna = True)
                cloudHeight = dataDraft["S_Ncbb"].mean(dim=["xt","yt", "time"], skipna = True)
                iceHeight  = dataDraft["S_Nibb"].mean(dim=["xt","yt", "time"], skipna = True)

                print("aeroHeight shape", aeroHeight.shape)

                aeroBin = aeroHeight[bini,:]
                cloudBin = cloudHeight[bini,:]
                iceBin = iceHeight[bini,:]

                totalBin = aeroBin + cloudBin + iceBin

                aeroFrac = aeroBin/totalBin
                cloudFrac = cloudBin/totalBin
                iceFrac = iceBin/totalBin



                aeroFrac.plot(ax=ax, color = aeroColor, y = "zt")
                cloudFrac.plot(ax=ax, color = cloudColor, y = "zt")
                iceFrac.plot(ax=ax, color = iceColor, y = "zt")


                bininame = str(bini +1)
                if True:
                    label = " ".join([drafType, "Bin", bininame])
                    facecolor = "black"

                    PlotTweak.setArtist(ax, {label:facecolor}, loc = (0.025, 0.85), framealpha = 0.8)

                    if axIndex == 0:
                        collectionOfLabelsColors = {"Aerosol": aeroColor, "Cloud": cloudColor, "Ice": iceColor}
                        PlotTweak.setArtist(ax, collectionOfLabelsColors, ncol = 3, loc = (0.75,1.12))

                ##############
                ax.set_title("")

                PlotTweak.setYLim(ax, end = yend)
                yticks = PlotTweak.setYticks(ax, end = yend, interval = 100)
                shownYLabelsBoolean = PlotTweak.setYLabels(ax, yticks, end = yend, interval = 200)
                PlotTweak.setYTickSizes(ax, shownYLabelsBoolean)


                PlotTweak.setXLim(ax,  end = xaxisend)
                xticks = PlotTweak.setXticks(ax, end =xaxisend, interval = 0.1, integer = False)
                xShownLabelsBoolean = PlotTweak.setXLabels(ax, xticks, end = xaxisend, interval = 0.2, integer=False)
                PlotTweak.setXTickSizes(ax, xShownLabelsBoolean)


                PlotTweak.setYaxisLabel(ax,"")
                if (bini == packing -1):
                    PlotTweak.setXaxisLabel(ax,"Time", "h")
                else:
                    PlotTweak.setXaxisLabel(ax,"")
                    PlotTweak.hideXTickLabels(ax)

                if draftIndex == 0:
                    PlotTweak.setXaxisLabel(ax,"Height", "m")

                if draftIndex == 1:
                    PlotTweak.hideYTickLabels(ax)

                if axIndex == 0:
                    ax.text( 0.5*xaxisend, yticks[-1]+yticks[1]*0.25, "Mean profile from " + PlotTweak.getUnitLabel(f"t_0={self.xstart}", "h")+ " to " +  PlotTweak.getUnitLabel(f"t_1={self.xend}", "h") + " limit: " + f"{self.draftLimit:.0e}"  , size=8)

            # end bini for loop
        # end draftIndex for loop
        fig.save()




def main(folder = os.environ["SIMULATIONFIGUREFOLDER"], datafolder = "/home/aholaj/Data/BinnedData", xstart = 2.1, xend = 33.0, draftLimit = 1e-3):

    figObject = ManuscriptFigures(folder, datafolder, xstart, xend, draftLimit)


    if False:
        figObject.figureUpdraftTimeseries()
    if True:
        figObject.figureUpdraftProfile()

if __name__ == "__main__":
    start = time.time()
    try:
        folder = sys.argv[1]
        datafolder = sys.argv[2]
        xstart = float(sys.argv[3])
        xend= float(sys.argv[4])
        draftLimit = float(sys.argv[5])
    except IndexError:
        folder = os.environ["SIMULATIONFIGUREFOLDER"]
        datafolder = "/home/aholaj/Data/BinnedData"
        xstart = 2.1
        xend = 33.0
        draftLimit = 1e-2

    main(folder =folder, datafolder = datafolder, xstart = xstart, xend = xend, draftLimit = draftLimit)

    end = time.time()
    print("Script completed in " + str(round((end - start),0)) + " seconds")
