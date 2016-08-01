
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

import gld.infra.Node.NodeStatistics;
import gld.infra.RoaduserFactory;

import java.awt.Canvas;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.Color;
import java.awt.event.*;

/**
 *
 * An <code>ExtendedTrackingView</code> shows a certain tracking graph.
 * To be used in combination with ExtendedTrackingController. Makes it 
 * possible to draw a graph for each concrete roaduser type and all roadusers.
 *
 * @author  Group GUI
 * @version 1.0
 */

public abstract class ExtendedTrackingView extends TrackingView
{
	protected boolean allTime;
	
  /**
   * Creates a new <code>ExtendedTrackingView</code>. A <code>TrackingView</code> is a component
   * displaying a tracking graph.
   *
   * @param _startCycle The cycle this view starts tracking.
   */
  public ExtendedTrackingView(int _startCycle)
  {
  	super(_startCycle);
  }
  

	/*============================================*/
	/* Basic GET and SET methods                  */
	/*============================================*/

	public void setAllTime(boolean b) { allTime = b; }
	public boolean getAllTime() { return allTime; }
	
  /** Returns the colors to be used when drawing the tracking graphs. */
  protected Color[] getColors() {
  	Color[] ruColors = new Color[MAX_TRACK];
  	ruColors[0] = Color.black;
  	for(int i=1; i<MAX_TRACK; i++)
  		ruColors[i] = RoaduserFactory.getColorByType(RoaduserFactory.statIndexToRuType(i));
  	return ruColors;
	}
	
	protected int getMaxTrack() { return RoaduserFactory.statArrayLength(); }
	
	protected String getSourceDesc(int i) { return RoaduserFactory.getDescByStatIndex(i); }

	public boolean useModes() { return true; }
}