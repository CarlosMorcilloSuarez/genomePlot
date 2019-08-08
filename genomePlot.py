#!/usr/bin/env python3
# -*- coding: utf-8 -*-

''' genomePlot.py	
        Library to create plots of genome regions
'''

__author__ = "Carlos Morcillo-Suarez"
__license__ = "GPL"
__version__ = "2019/08/08 12:03" # YYYY/MM/DD HH:MM
__email__ = "carlos.morcillo.upf.edu@gmail.com"
__copyright__ = "Copyright 2019, Carlos Morcillo-Suarez"

import matplotlib.pyplot as plt

def defineXTicks(initPosition, endPosition, modifyUnit = False):
    ''' 
        Defines the xTick values appropiate for the genome range of
        the plot
    '''
    segmentLength = endPosition - initPosition
    
    offset = 0
    for multiplicator in [1,10,100,1000,10000,100000,1000000,10000000,100000000]:
        for base in [1,2,5]:
            if segmentLength / (multiplicator * base) < 8:
                offset =  multiplicator * base
                firstTick = int(((initPosition/offset)+1))*offset
                lastTick = int((endPosition/offset))*offset
                ticks = range(firstTick,lastTick+1,offset)

                # Calculates unit
                if modifyUnit:
                    if 0 < offset/1000000:
                        unit = " (Mb)"
                        labels = [tick/1000000 for tick in ticks]
                    elif 0 < offset/1000:
                        unit = " (Kb)"
                        labels = [tick/1000 for tick in ticks]
                else:
                    unit = ""
                    labels = [tick for tick in ticks]
                    
                # Format labels
                labels = ["{:,}".format(label) for label in labels]
                
                return(ticks,labels, unit)
    

def createsGenomePlot(   ax,
                        chromosome = "", plotID = "", 
                        xOrigin = "", xWidth = "", 
                        yScale = 0, yLabel = ""
                        ):
                        
        # Scales and deletes frames
        ax.axis([xOrigin,xOrigin+xWidth,0-int(yScale/20),yScale])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        

        # Title
        ax.set_title(plotID + "  CHR: " + chromosome,
                pad = 20,
                #y = 1.05,
                x = -0.05,
                ha = 'left',
                color = "black",
                style = "italic"
                )

        # y axis
        ax.set_ylabel(yLabel+"\n")
                
        # x axis
        xtickPositions, labels, unit = defineXTicks(xOrigin,xOrigin+xWidth)
        
        ax.set_xlabel('Position'+unit)
                            
        ax.set_xticks(xtickPositions)
        ax.set_xticklabels(labels) 
            
    
def createGenomePlotWithPoints( ax, positions, values, 
                                chromosome = "", plotID = "", 
                                xOrigin = "", xWidth = "", 
                                yScale = 0, yLabel = "",
                                color = "blue",
                                markersize = 3):    
    # Defines x and y scale of plot    
    if xWidth == "":
        xWidth = max(positions)-min(positions)
    xWidth = int(xWidth*1.01)
    if xOrigin == "": 
        xOrigin = min(positions)
    xOrigin = int(xOrigin - (xWidth*0.005))
    if yScale == 0:
        yScale = int(max(values)*1.05)
        
    # Create Genome Plot
    createsGenomePlot(   ax,
                        chromosome, plotID, 
                        xOrigin, xWidth, 
                        yScale, yLabel,
                        )
        
    # plots values
    addPointsToGenomePlot(ax, positions,values,color,markersize)


def createGenomePlotWithSegments( ax, startPositions, endPositions, values, 
                                chromosome = "", plotID = "", 
                                xOrigin = "", xWidth = "", 
                                yScale = 0, yLabel = "",
                                color = "blue"):    
    # Defines x and y scale of plot    
    if xWidth == "":
        xWidth = max(endPositions)-min(startPositions)
    xWidth = int(xWidth*1.01)
    if xOrigin == "": 
        xOrigin = min(startPositions)
    xOrigin = int(xOrigin - (xWidth*0.005))
    if yScale == 0:
        yScale = int(max(values)*1.05)
        
    # Create Genome Plot
    createsGenomePlot(  ax,
                        chromosome, plotID, 
                        xOrigin, xWidth, 
                        yScale, yLabel,
                        )
        
    # plots values
    addSegmentsToGenomePlot(ax, startPositions, endPositions, values,color)     


def addPointsToGenomePlot(ax, positions,values,color,markersize = 3):                    
    # plots values
    ax.plot(positions,
            values,
            "o",
            color = color,
            markersize = markersize)

    # plots depth values that are higher that yScale
    # as red dots with yScale value
    yScale =  ax.get_ylim()[1]
    values = [None if value == None or value < yScale
                    else yScale
                    for value
                    in values
                    ]    
    ax.plot(positions,
            values,
            "ro",
            markersize = 8) 
    
    
def addSegmentsToGenomePlot(ax, startPositions, endPositions, values, color,
                            linewidth = 3):
    # plots values
    data = []
    for startPosition, endPosition, value in zip(startPositions, 
                                                    endPositions, 
                                                    values):
        # Corrects values too high
        yScale =  ax.get_ylim()[1]
        if yScale < value:
            value = yScale*0.99
            currentColor = "red"
        else:
            currentColor = color
            
        data.append((startPosition, endPosition))
        data.append((value,value))
        data.append(currentColor)
    
    ax.plot(*data,linewidth=linewidth)
    
    
def mrCaNaVaRPlot(ax, startPositions, endPositions, values, 
                    chromosome, plotID, xOrigin, xWidth, yScale):
    '''
        mrCaNaNaRPlot takes values according to *bed files format:
            first position of Genome = 0
            segment spands from initPosition to endPosition-1
    '''
    # Change to 0-based to 1-based genome annotation
    startPositions = [startPosition+1 for startPosition in startPositions]
    # Plots Segments
    createGenomePlotWithSegments( ax, startPositions, endPositions, values, 
                                chromosome = chromosome, 
                                plotID = plotID, 
                                xOrigin = xOrigin,
                                xWidth = xWidth,
                                yScale = yScale,
                                yLabel = "mrCaNaVaR Calling"
                                )
    
    ax.spines['bottom'].set_visible(False)
    
    # Adds horizontal lines
    ax.axhline(y=0,linewidth=1, color='gray')
    ax.axhline(y=2,linewidth=1, color='gray')
    ax.axhline(y=4,linewidth=1, color='gray')
    ax.axhline(y=10,linewidth=1, color='gray')
    
