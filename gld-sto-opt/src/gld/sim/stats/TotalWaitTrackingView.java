
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

/**
 *
 * TrackingView that tracks the average trip waiting time.
 *
 * @author  Group GUI
 * @version 1.0
 */

public class TotalWaitTrackingView extends ExtendedTrackingView
{
	NodeStatistics[][] stats;
	
	public TotalWaitTrackingView(int startCycle, SimModel model)
	{
		super(startCycle);
		stats = model.getInfrastructure().getEdgeNodeStatistics();
	}

	/** Returns the next sample to be 'tracked'. */
	protected float nextSample(int index) 
	{ 
		float sample = 0, count = 0;
		int ru;

		for(int i=0; i<stats.length; i++) {
			if(allTime)
				ru = stats[i][index].getTotalRoadusers();
			else
				ru = Math.min(Node.STAT_NUM_DATA, stats[i][index].getTotalRoadusers());
			sample += stats[i][index].getAvgWaitingTime(allTime) * ru;
			count += ru;
		}

		return count == 0 ? 0 : sample / count;
	}

	protected String getYLabel() { return "delay (cycles)"; }

	public String getDescription() { return "average trip waiting time"; }
}