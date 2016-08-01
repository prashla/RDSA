
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

import gld.infra.Junction;
import gld.infra.Node.NodeStatistics;

/**
 *
 * TrackingView that tracks the waiting queue length of one EdgeNode.
 *
 * @author  Group GUI
 * @version 1.0
 */

public class JunctionWaitTrackingView extends ExtendedTrackingView
{
	NodeStatistics[] stats;
	int id;
	
  public JunctionWaitTrackingView(int startCycle, Junction junction)
  {
		super(startCycle);
		stats = junction.getStatistics();
		id = junction.getId();
  }

	/** Returns the next sample to be 'tracked'. */
	protected float nextSample(int src) 
	{ 
		return stats[src].getAvgWaitingTime(allTime);
	}
	
	/** Returns the description for this tracking window. */
	public String getDescription() { return "junction " + id + " - average junction waiting time"; }
	
	protected String getYLabel() { return "delay (cycles)"; }
}