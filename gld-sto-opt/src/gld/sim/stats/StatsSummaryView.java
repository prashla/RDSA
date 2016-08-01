
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

/**
*
* Extension of StatisticsView showing a summary of all data.
*
* @author Group GUI
* @version 1.0
*/

public class StatsSummaryView extends StatisticsView
{
	protected final static int LINE_HEIGHT = 20;
	protected int x, y;
	
	public StatsSummaryView(StatisticsController parent, StatisticsModel stats)
	{
		super(parent, stats);
		parent.setScrollMax(0,0);
	}
	
	public void paintStats(Graphics g)
	{
		g.setColor(Color.black);
		g.setFont(infoFont);

		y = paintArea.y + LINE_HEIGHT;
				
		infoLine(g, "Nodes: " + stats.getNumSpecialNodes() + " special nodes, " + stats.getNumJunctions() + " junctions");
		emptyLine();
		infoLine(g, "Total number of roadusers that has arrived at its destination: " + stats.getRoadusersArrived());
		infoLine(g, "Total number of junction crossings: " + stats.getJunctionCrossings());
		emptyLine();
		infoLine(g, "Average trip waiting time (based on all roadusers arrived): " + stats.getAllTimeTripWT());
		if(stats.getLastXTripCount() != stats.getRoadusersArrived())
				infoLine(g, "Average trip waiting time (based on last "+stats.getLastXTripCount()+" roadusers arrived): " + stats.getLastXTripWT());
		infoLine(g, "Average junction waiting time (based on all junction crossings): " + stats.getAllTimeJunctionWT());
		if(stats.getLastXJunctionCount() != stats.getJunctionCrossings())
				infoLine(g, "Average junction waiting time (based on last "+stats.getLastXJunctionCount()+" junction crossings): " + stats.getLastXJunctionWT());
	}
	
	protected void infoLine(Graphics g, String s) {
		g.drawString(s, x, y);
		y += LINE_HEIGHT;
	}
	
	protected void emptyLine() { y += LINE_HEIGHT; }

	protected void paintAreaChanged() {
		x = paintArea.x;
	}
}