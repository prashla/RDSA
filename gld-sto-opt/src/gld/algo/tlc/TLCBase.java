/*-----------------------------------------------------------------------
 * Copyright (C) 2001 Green Light District Team, Utrecht University 
 *
 * This program (Green Light District) is free software.
 * You may redistribute it and/or modify it under the terms
 * of the GNU General Public License as published by
 * the Free Software Foundation (version 2 or later).
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * See the documentation of Green Light District for further information.
 *------------------------------------------------------------------------*/

package gld.algo.tlc;

import gld.*;
import gld.sim.*;
import gld.algo.tlc.*;
import gld.infra.*;
import gld.utils.*;
import gld.xml.*;

import java.awt.*;
import java.awt.event.*;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Dictionary;
import java.util.List;
import java.util.Random;

/**
 * 
 * This is the abstract class for Traffic light algorithms. It is informed about
 * every movement made by road users. In this way not every road user has to be
 * iterated. By using this information it provides a table containing
 * Q-values(reward values) for each trafficlight in it's 'Green' setting.
 * 
 * @author Group Algorithms
 * @version 1.0
 */
public abstract class TLCBase extends TLController
{
    // For prioritizing main road lanes
    File lanesFile;
    int priorityLanes[];
    double laneWeights[];

    ListCombinationGenerator<Integer> actionList;

    int[] queueLengths;
    int[] elapsedTimes;

    public static final int RED = 0;
    public static final int GREEN = 1;
    
    public static final int MAX_RED_CYCLES = 300;

    
    // In case you want to run a RL algo till it converges and then fx the policy, use this constant
    public static final int RUNLENGTH = 3000;
    Random randGen;

    // TLC vars
    protected Infrastructure infrastructure;
    protected TrafficLight[][] tls;
    protected Node[] nodes;
    protected Node[] allnodes;
    protected Junction[] junctions;
    protected int num_nodes;
    protected int numRealSigns = 0;
    protected double cost;
 
    /**
     * The constructor for TL controllers
     * 
     * @param The
     *            infrastructure being used.
     */
    TLCBase()
    {
    }

    TLCBase(Infrastructure i)
    {
        setInfrastructure(i);
        infra = i;
        randGen = new Random(1015);
    }

    
    public double getAlphaWeight() { return alphaWeight; }
    public void setAlphaWeight(double k) { alphaWeight = k;}    


    // Get the lanes with higher priority from lanes.txt
    public void setLanesFile(String k)
    {
        lanesFile = new File(k);
        String fileContents = Fileutils.getContents(lanesFile);
        //System.out.println("Original file contents: " + fileContents);
        String x[] = fileContents.split(" ");
        priorityLanes = new int[x.length];
        // Get lane ids
        nodes = infra.getAllNodes();
        int numRealSigns = 0;
        for (int i = 0; i < nodes.length; i++)
        {
            numRealSigns += nodes[i].getNumSigns();
        }

        laneNumbers = new int[numRealSigns];
        laneWeights = new double[numRealSigns];
        int sPos = 0;
        int num_lanes;
        int num_nodes = nodes.length;;
        for (int i = 0; i < num_nodes; i++)
        {
            num_lanes = tld[i].length;
            for (int j = 0; j < num_lanes; j++)
            {
                if (tld[i][j].getTL().getType() == Sign.TRAFFICLIGHT)
                    laneNumbers[sPos++] = tld[i][j].getTL().getLane().getId();
            }
        }

        for (int i = 0; i < x.length; i++)
        {
            priorityLanes[i] = Integer.parseInt(x[i]);
        }

        for (int i = 0; i < laneNumbers.length; i++)
        {
            boolean found = false;
            for (int j = 0; j < priorityLanes.length; j++)
            {
                if (priorityLanes[j] == laneNumbers[i])
                {
                    found = true;
                    break;
                }
            }
            if (found)
                laneWeights[i] = alphaWeight;
            else
                laneWeights[i] = 1 - alphaWeight;
        }
    }

    public void initFiles(int simNr)
    {
    }
    
    public void updateNumElapsedCycles()
    {
        TLDecision curDec;
        // System.out.println("NumGreenCycles array:");
        for (int i = 0; i < tld.length; i++)
        { // for all nodes
            // System.out.println("node " + i);
            for (int j = 0; j < tld[i].length; j++)
            { // for all inbound lanes
                // in node

                curDec = tld[i][j];
                if (tld[i][j].getTL().getState())
                    curDec.addNumGreenCycles(1);
                if (!tld[i][j].getTL().getState())
                    curDec.addNumRedCycles(1);
                if(curDec.getNumRedCycles() > MAX_RED_CYCLES) {curDec.setNumRedCycles(MAX_RED_CYCLES);} 
                // System.out.print(curDec.getNumGreenCycles() + " ");
            }
        }
    }

    public void clearNumGreenCycles()
    {
        for (int i = 0; i < tld.length; i++)
        { // for all nodes
            for (int j = 0; j < tld[i].length; j++)
            { // for all inbound lanes
                // in node
                if (tld[i][j].getTL().getCycleSwitchedToRed() == getCurCycle() - 1)
                {
                    tld[i][j].setNumGreenCycles(0);
                }
            }
        }
    }

    public void clearNumRedCycles()
    {
        for (int i = 0; i < tld.length; i++)
        { // for all nodes
            for (int j = 0; j < tld[i].length; j++)
            { // for all inbound lanes
                // in node
                if (tld[i][j].getTL().getCycleSwitchedToGreen() == getCurCycle() - 1)
                {
                    tld[i][j].setNumRedCycles(0);
                }
            }
        }
    }

}