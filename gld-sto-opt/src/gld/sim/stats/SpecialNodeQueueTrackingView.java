
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

/**
 *
 * TrackingView that tracks the waiting queue length of one EdgeNode.
 *
 * @author  Group GUI
 * @version 1.0
 */

public class SpecialNodeQueueTrackingView extends ExtendedTrackingView
{	SpecialNode node;
	
  public SpecialNodeQueueTrackingView(int startCycle, SpecialNode node)
  {	super(startCycle);
		this.node=node;
  }

	/** Returns the next sample to be 'tracked'. */
	protected float nextSample(int src) 
	{ 	return node.getWaitingQueueLength();
	}
	
	/** Returns the description for this tracking window. */
	public String getDescription() { return "special node " + node.getId() + " - waiting queue length"; }
	
	protected String getSourceDesc(int src) { return "length"; }
	protected String getYLabel() { return "queue length (roadusers)"; }
}
