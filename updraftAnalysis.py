#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 10:29:56 2020

@author: Jaakko Ahola, Finnish Meteorological Institute
@licence: MIT licence Copyright
"""

import copy
import matplotlib
import numpy
import time
import sys
import os
import pathlib
import xarray

sys.path.append("../LES-03plotting")
from Figure import Figure
from PlotTweak import PlotTweak
from Colorful import Colorful
from SimulationDataAnalysis import SimulationDataAnalysis
from Simulation import Simulation
matplotlib.use("PDF")

class ManuscriptFigures:

    def __init__(self, figurefolder, datafolder, xstart = 28., xend = 32., draftLimit = 0.1):

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
        self.simulation.setAUXDataset("Updraft", "case_isdac_LVL5_3D_iceD_inter_48h_w.nc")


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
                    draftType = "Up-draft"
                else:
                    dataDraft = datasetHeight.where(datasetHeight["w"] < - self.draftLimit, drop=True)
                    draftType = "Down-draft"

                for bini in range(packing):
                    print("")
                    print("height", realHeight, draftType,"bini", bini)
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

                    bininame = str(bini +1)
                    if True:
                        label = " ".join([draftType, "Bin", bininame + ",", "Total", r"$N_0$",  str(int(pointZero)) + ",", "\nMin", r"$N$", str(int(pointEnd)), "$(kg^{-1})$"  ])
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
        if xarray.__version__ == "0.14.1":
            updraft.yt.values = updraft.yt.values-25.
        else:
            updraft = updraft.assign_coords(yt = (updraft.yt.values-25.))

        dataset = xarray.merge([aero,cloud,ice,updraft])

        aeroColor = Colorful.getDistinctColorList("red")
        cloudColor = Colorful.getDistinctColorList("navy")
        iceColor = Colorful.getDistinctColorList("cyan")
        totalColor = Colorful.getDistinctColorList("green")

        xaxisend = 1.0
        ystart = 400.
        yend = 850.

        loopedBins = range(packing)[1:3]

        fig = Figure(self.figurefolder,"figure8",
                    figsize = [4.724409448818897, 3],
                     ncols = 2, nrows =len(loopedBins),
                     wspace=0.06, left=0.12, bottom = 0.12, top=0.85, right=0.98)
        print("figsize", fig.getFigSize())

        timeBegin = dataset.time.sel(time=self.xstart, method = "nearest").item()
        timeEnd = dataset.time.sel(time=self.xend, method = "nearest").item()
        datasetTime = dataset.sel(time = slice(timeBegin,timeEnd))


        subFigureName =  [chr(ord('a') + i) for i in range(packing*2)]

        for draftIndex in range(2):
            if draftIndex == 0:
                dataDraft = datasetTime.where(datasetTime["w"] > self.draftLimit, drop=True)
                draftType = "Up-draft"
            else:
                dataDraft = datasetTime.where(datasetTime["w"] < - self.draftLimit, drop=True)
                draftType = "Down-draft"

            biniCounter = 0
            for bini in loopedBins:
                print("")
                print(draftType, "timeBegin", timeBegin, "timeEnd", timeEnd)
                axIndex =  biniCounter*2 + draftIndex

                ax = fig.getAxes(axIndex)

                if dataDraft["w"].size == 0:
                    ax.axis("off")
                    continue

                aeroHeight = dataDraft["S_Nabb"].mean(dim=["xt","yt", "time"], skipna = True)
                cloudHeight = dataDraft["S_Ncbb"].mean(dim=["xt","yt", "time"], skipna = True)
                iceHeight  = dataDraft["S_Nibb"].mean(dim=["xt","yt", "time"], skipna = True)


                aeroBin = aeroHeight[bini,:]
                cloudBin = cloudHeight[bini,:]
                iceBin = iceHeight[bini,:]

                totalBin = aeroBin + cloudBin + iceBin

                aeroFrac = aeroBin/totalBin
                cloudFrac = cloudBin/totalBin
                iceFrac = iceBin/totalBin

                pointMax = totalBin.max()


                totalBinRelative  = totalBin / pointMax.item()
                totalBinRelative.plot(ax = ax, color = totalColor, y = "zt")


                aeroFrac.plot(ax=ax, color = aeroColor, y = "zt")
                cloudFrac.plot(ax=ax, color = cloudColor, y = "zt")
                iceFrac.plot(ax=ax, color = iceColor, y = "zt")

                if packing < 7:
                    if bini == (packing - 1):
                        bininame = str(bini + 1 ) + " - 7"
                    else:
                        bininame = str(bini +1)
                else:
                    bininame = str(bini +1)

                ##############
                ax.set_title("")

                PlotTweak.setYLim(ax, start = ystart, end = yend)
                yticks = PlotTweak.setYticks(ax, start = ystart, end = yend, interval = 50)
                shownYLabelsBoolean = PlotTweak.setYLabels(ax, yticks, end = yend, interval = 100, integer = False)
                PlotTweak.setYTickSizes(ax, shownYLabelsBoolean)


                PlotTweak.setXLim(ax,  end = xaxisend)
                xticks = PlotTweak.setXticks(ax, end =xaxisend, interval = 0.1, integer = False)
                xTickLabels = copy.deepcopy(xticks)
                for xtickInd, xtickValue in enumerate(xTickLabels):
                    if (xtickInd == 0) or (xtickValue == xTickLabels[-1]):
                        xTickLabels[xtickInd] = f"{xtickValue:.0f}"
                    else:
                        xTickLabels[xtickInd] = f"{xtickValue:.1f}"

                xShownLabelsBoolean = PlotTweak.setXLabels(ax, xticks, end = xaxisend, interval = 0.2, integer=False)
                ax.set_xticklabels(xTickLabels)
                PlotTweak.setXTickSizes(ax, xShownLabelsBoolean)

                xlimits = ax.get_xlim()
                ylimits = ax.get_ylim()
                PlotTweak.setXaxisLabel(ax,"")
                PlotTweak.setYaxisLabel(ax,"")

                print("bini", biniCounter, "loopedBins -1", len(loopedBins)-1)
                if (biniCounter != len(loopedBins) -1):
                     PlotTweak.hideXTickLabels(ax)


                if draftIndex == 0:
                    PlotTweak.setYaxisLabel(ax,"Height", "m")

                if draftIndex == 1:
                    PlotTweak.hideYTickLabels(ax)

                if axIndex == 0:
                    ax.text( 0.05*(xlimits[1]-xlimits[0]), (ylimits[1]-ylimits[0])*0.05+ylimits[1],
                            "Mean profile from " + PlotTweak.getUnitLabel(f"t_0={self.xstart}", "h")+\
                                " to " +  PlotTweak.getUnitLabel(f"t_1={self.xend}", "h") +\
                                    " limit: " + PlotTweak.getUnitLabel(f"{self.draftLimit:.1f}", "m\ s^{-1}")
                            , size=8)
                if True:
                    label = " ".join([subFigureName[axIndex] + ")", draftType, "Bin", bininame])

                    PlotTweak.setAnnotation(ax, label, xPosition = 0.05*(xlimits[1]-xlimits[0])+xlimits[0], yPosition = 0.1*(ylimits[1]-ylimits[0])+ylimits[0])

                    if axIndex == 0:
                        collectionOfLabelsColors = {"Aerosol": aeroColor, "Cloud": cloudColor, "Ice": iceColor, "Total":totalColor}
                        PlotTweak.setArtist(ax, collectionOfLabelsColors, ncol = 4, loc = (0.10,1.19))
                biniCounter += 1

            # end bini for loop
        # end draftIndex for loop
        fig.save()


    def figureUpdraftProfileLog(self):
        packing = 4


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

        ystart = 400.
        yend = 850.

        fig = Figure(self.figurefolder,"figureProfileDraftsLog",
                    figsize = [4.724409448818897, 2.5],
                     ncols = 2, nrows =1,
                     wspace=0.06, left=0.12, bottom = 0.12, top=0.86, right=0.98)
        print("figsize", fig.getFigSize())

        timeBegin = dataset.time.sel(time=self.xstart, method = "nearest").item()
        timeEnd = dataset.time.sel(time=self.xend, method = "nearest").item()
        datasetTime = dataset.sel(time = slice(timeBegin,timeEnd))


        subFigureName =  [chr(ord('a') + i) for i in range(packing*2)]
        loopedBins = range(packing)[1:2]
        for draftIndex in range(2):
            if draftIndex == 0:
                dataDraft = datasetTime.where(datasetTime["w"] > self.draftLimit, drop=True)
                draftType = "Up-draft"
            else:
                dataDraft = datasetTime.where(datasetTime["w"] < - self.draftLimit, drop=True)
                draftType = "Down-draft"

            biniCounter = 0
            for bini in loopedBins:
                print("")
                print(draftType, "timeBegin", timeBegin, "timeEnd", timeEnd)
                axIndex =  biniCounter*2 + draftIndex
                biniCounter += 1
                ax = fig.getAxes(axIndex)

                if dataDraft["w"].size == 0:
                    ax.axis("off")
                    continue

                aeroHeight = dataDraft["S_Nabb"].mean(dim=["xt","yt", "time"], skipna = True)
                cloudHeight = dataDraft["S_Ncbb"].mean(dim=["xt","yt", "time"], skipna = True)
                iceHeight  = dataDraft["S_Nibb"].mean(dim=["xt","yt", "time"], skipna = True)


                aeroBin = aeroHeight[bini,:]
                cloudBin = cloudHeight[bini,:]
                iceBin = iceHeight[bini,:]

                aeroFrac = numpy.log10(aeroBin)
                cloudFrac = numpy.log10(cloudBin)
                iceFrac = numpy.log10(iceBin)


                aeroFrac.plot(ax=ax, color = aeroColor, y = "zt")
                cloudFrac.plot(ax=ax, color = cloudColor, y = "zt")
                iceFrac.plot(ax=ax, color = iceColor, y = "zt")

                if packing < 7:
                    if bini == (packing - 1):
                        bininame = str(bini + 1 ) + " - 7"
                    else:
                        bininame = str(bini +1)
                else:
                    bininame = str(bini +1)

                ##############
                ax.set_title("")

                PlotTweak.setYLim(ax, start = ystart, end = yend)
                yticks = PlotTweak.setYticks(ax, start = ystart, end = yend, interval = 50)
                shownYLabelsBoolean = PlotTweak.setYLabels(ax, yticks, end = yend, interval = 100, integer = False)
                PlotTweak.setYTickSizes(ax, shownYLabelsBoolean)


                xlimits = ax.get_xlim()
                ylimits = ax.get_ylim()
                print("limits", xlimits, ylimits)

                PlotTweak.setXaxisLabel(ax,"")
                PlotTweak.setYaxisLabel(ax,"")

                if draftIndex == 0:
                    PlotTweak.setYaxisLabel(ax,"Height", "m")

                if draftIndex == 1:
                    PlotTweak.hideYTickLabels(ax)

                if axIndex == 0:
                    ax.text( 0.1*(xlimits[1]-xlimits[0])+xlimits[0], (ylimits[1]-ylimits[0])*0.05+ylimits[1],
                            "Mean profile from " + PlotTweak.getUnitLabel(f"t_0={self.xstart}", "h")+\
                                " to " +  PlotTweak.getUnitLabel(f"t_1={self.xend}", "h") +\
                                    " limit: " + f"{self.draftLimit:.0e}"  , size=8)
                if True:
                    label = " ".join([subFigureName[axIndex] + ")", draftType, "Bin", bininame])

                    PlotTweak.setAnnotation(ax, label,
                                            xPosition = 0.08*(xlimits[1]-xlimits[0])+xlimits[0],
                                            yPosition = 0.08*(ylimits[1]-ylimits[0])+ylimits[0])

                    if axIndex == 0:
                        collectionOfLabelsColors = {"Aerosol": aeroColor, "Cloud": cloudColor, "Ice": iceColor}
                        PlotTweak.setArtist(ax, collectionOfLabelsColors, ncol = 3, loc = (0.15,1.19))

            # end bini for loop
        # end draftIndex for loop
        fig.save()



def main(folder = os.environ["SIMULATIONFIGUREFOLDER"], datafolder = "/home/aholaj/Data/BinnedData", xstart = 28., xend = 32., draftLimit = 0.1):

    figObject = ManuscriptFigures(folder, datafolder, xstart, xend, draftLimit)


    if False:
        figObject.figureUpdraftTimeseries()
    if True:
        figObject.figureUpdraftProfile()
    if False:
        figObject.figureUpdraftProfileLog()


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
        datafolder = pathlib.Path( os.environ["SIMULATIONDATAROOTFOLDER"] ) / "case_isdac_LVL5_3D_iceD_inter_48h"
        xstart = 28.0
        xend = 32.0
        draftLimit = 0.1

    main(folder =folder, datafolder = datafolder, xstart = xstart, xend = xend, draftLimit = draftLimit)

    end = time.time()
    print("Script completed in " + str(round((end - start),0)) + " seconds")
