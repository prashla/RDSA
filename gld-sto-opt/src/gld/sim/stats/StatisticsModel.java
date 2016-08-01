
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

package gld.sim.stats;

import gld.infra.*;
import gld.infra.Node.NodeStatistics;
import gld.sim.SimModel;

import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.io.*;

/**
*
* Class to maintain statistics as shown in extensions of StatisticsView.
*
* @author Group GUI
* @version 1.0
*/

public class StatisticsModel extends Observable
{
	/** Separator used when saving data. */
	public static String SEP = "\t";
	/** The length of the delay tables maintained in each Node. */
	protected final static int statNumData = Node.STAT_NUM_DATA;
	/** Determines whether to use all time averages or 'last <statNumData> averages' */
	protected boolean allTimeAvg;
	/** The SimModel to get the statistical data from. */
	protected SimModel model;

	// All variables used to maintain the statistics.
	protected String infraName, infraAuthor, simName;
	protected int numNodes, numSpecial, numJunctions, cycle;
	protected NodeStatistics[][] nodeStats;
	protected float allTimeTripWT, lastXTripWT, allTimeJunctionWT, lastXJunctionWT;
	protected int roadusersArrived, junctionCrossings, lastXTripCount, lastXJunctionCount;


	/**
	* Creates a <code>StatisticsModel</code>.
	*
	* @param _model The <code>SimModel</code> statistics should be read from.
	*/
	public StatisticsModel(SimModel _model)
	{
		model = _model;
		refresh();
	}





	/*============================================*/
	/* GET AND SET                                */
	/*============================================*/

	/** Sets the SimModel to be used. */
	public void setSimModel(SimModel _model) { model = _model; refresh(); }
	/** Returns the SimModel to be shown. */
	public SimModel getSimModel() { return model; }
	/** Sets whether all time averages should be used or not. */
	public void setAllTimeAvg(boolean b) { allTimeAvg = b; setChanged(); notifyObservers(); }
	/** Returns whether all time averages should be used or not. */
	public boolean getAllTimeAvg() { return allTimeAvg; }

	/// GET ///
	public String getInfraName() { return infraName; }
	public String getInfraAuthor() { return infraAuthor; }
	public String getSimName() { return simName; }
	public int getNumNodes() { return numNodes; }
	public int getNumSpecialNodes() { return numSpecial; }
	public int getNumJunctions() { return numJunctions; }
	public int getCycle() { return cycle; }
	public NodeStatistics[][] getNodeStatistics() { return nodeStats; }
	public float getAllTimeTripWT() { return allTimeTripWT; }
	public float getLastXTripWT() { return lastXTripWT; }	
	public float getAllTimeJunctionWT() { return allTimeJunctionWT; }
	public float getLastXJunctionWT() { return lastXJunctionWT; }
	public int getRoadusersArrived() { return roadusersArrived; }
	public int getJunctionCrossings() { return junctionCrossings; }
	public int getLastXTripCount() { return lastXTripCount; }
	public int getLastXJunctionCount() { return lastXJunctionCount; }



	/*============================================*/
	/* REFRESHING                                 */
	/*============================================*/

	/** Refresh the statistical data from the model and repaint. */
	public void refresh()
	{
		Infrastructure infra = model.getInfrastructure();
		infraName = infra.getTitle();
		infraAuthor = infra.getAuthor();
		simName = model.getSimName();
		cycle = model.getCurCycle();
		
		numNodes = infra.getNumNodes();
		numSpecial = infra.getNumSpecialNodes();
		numJunctions = infra.getNumJunctions();
		NodeStatistics[][] ns = infra.getNodeStatistics();
		nodeStats = new NodeStatistics[ns.length][];
		for(int i=0; i<ns.length; i++) {
			nodeStats[i] = new NodeStatistics[ns[i].length];
			for(int j=0; j<ns[i].length; j++)
				nodeStats[i][j] = (NodeStatistics)ns[i][j].getClone();
		}

		roadusersArrived = junctionCrossings = lastXTripCount = lastXJunctionCount = 0;
		allTimeTripWT = lastXTripWT = allTimeJunctionWT = lastXJunctionWT = 0;

		for(int i=0; i<numSpecial; i++) {
			int ru = nodeStats[i][0].getTotalRoadusers();
			roadusersArrived += ru;
			allTimeTripWT += nodeStats[i][0].getAvgWaitingTime(true) * ru;
			int tmp = Math.min(ru, statNumData);
			lastXTripWT += nodeStats[i][0].getAvgWaitingTime(false) * tmp;
			lastXTripCount += tmp;
		}
		allTimeTripWT = roadusersArrived > 0 ? allTimeTripWT / roadusersArrived : 0;
		lastXTripWT = lastXTripCount > 0 ? lastXTripWT / lastXTripCount : 0;

		for(int i=numSpecial; i<numNodes; i++) {
			int ru = nodeStats[i][0].getTotalRoadusers();
			junctionCrossings += ru;
			allTimeJunctionWT += nodeStats[i][0].getAvgWaitingTime(true) * ru;
			int tmp = Math.min(ru, statNumData);
			lastXJunctionWT += nodeStats[i][0].getAvgWaitingTime(false) * tmp;
			lastXJunctionCount += tmp;
		}
		allTimeJunctionWT = junctionCrossings > 0 ? allTimeJunctionWT / junctionCrossings : 0;
		lastXJunctionWT = lastXJunctionCount > 0 ? lastXJunctionWT / lastXJunctionCount : 0;
		
		setChanged();
		notifyObservers();
	}
	
	
	

	
	/*============================================*/
	/* SAVING                                     */
	/*============================================*/
	
	/**
	 * Save data to a CSV file. 
	 */
	protected void saveData(String filename) throws IOException
	{
		int[] ruTypes = new int[nodeStats[0].length];
		for(int i=0; i<ruTypes.length; i++)
			ruTypes[i] = RoaduserFactory.statIndexToRuType(i);
		
		PrintWriter out=new PrintWriter(new FileWriter(new File(filename)));
		out.println("# Data exported by Green Light District"); out.println("#");
		Infrastructure infra = model.getInfrastructure();
		out.println("# Infrastructure: \"" + infra.getTitle() + "\" by " + infra.getAuthor());
		out.println("# Simulation: \"" + model.getSimName() + "\""); out.println("#");
		out.println("# Data at cycle: " + cycle);
		out.println("# #nodes = " + numNodes + ", #specialnodes = " + numSpecial + ", #junctions = " + numJunctions);
		out.println("#"); 
		out.println("# Total number of roadusers that has arrived at its destination: " + roadusersArrived); 
		out.println("# Total number of junction crossings: " + junctionCrossings); 
		out.println("# Average trip waiting time (based on all roadusers arrived): " + allTimeTripWT); 
		if(roadusersArrived != lastXTripCount) 
			out.println("# Average trip waiting time (based on last " + lastXTripCount + " roadusers arrived): " + lastXTripWT); 
		out.println("# Average junction waiting time (based on all junction crossings): " + allTimeJunctionWT); 
		if(junctionCrossings != lastXJunctionCount)
			out.println("# Average junction waiting time (based on last " + lastXJunctionCount + " junction crossings): " + lastXJunctionWT); 
		out.println("#"); 
		out.println("#"); out.println("# EdgeNodes"); 
		out.println("# Data format: " + "<id"+SEP+"ruType"+SEP+"roadusersArrived"+SEP+"avgTripWaitingTimeAllTime"+SEP+"avgTripWaitingTimeLast" + Node.STAT_NUM_DATA + ">"); out.println("#");
		for(int id=0; id<numNodes; id++) {
			if(id==numSpecial) { out.println("#"); out.println("# Junctions"); out.println("#"); }
			for(int statIndex=0; statIndex<nodeStats[id].length; statIndex++) {
				NodeStatistics ns = nodeStats[id][statIndex];
				out.println(id + SEP + ruTypes[statIndex] + SEP + ns.getTotalRoadusers() + SEP + 
					ns.getAvgWaitingTime(true) + SEP + ns.getAvgWaitingTime(false));
			}
		}
		out.close();
	}	
}