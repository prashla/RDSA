
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

import gld.infra.SpecialNode;
import gld.sim.SimModel;

/**
 *
 * TrackingView that tracks the sum of the lengths of all waiting queues.
 *
 * @author  Group GUI
 * @version 1.0
 */

public class AllQueuesTrackingView extends ExtendedTrackingView
{
	SpecialNode[] specialNodes;
	int numSpecialNodes;
	
  public AllQueuesTrackingView(int startCycle, SimModel model)
  {
		super(startCycle);
		specialNodes = model.getInfrastructure().getSpecialNodes();
		numSpecialNodes = specialNodes.length;
  }

	/** Returns the next sample to be 'tracked'. */
	protected float nextSample(int src) 
	{ 
		int sample = 0;
		for(int i=0; i<numSpecialNodes; i++)
			sample += specialNodes[i].getWaitingQueueLength();
		return sample; 
	}
	
	/** Returns the description for this tracking window. */
	public String getDescription() { return "total waiting queue length"; }
	
	protected String getSourceDesc(int i) { return "length"; }
	protected String getYLabel() { return "queue length (roadusers)"; }
}
