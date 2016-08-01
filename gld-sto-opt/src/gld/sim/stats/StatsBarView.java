
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
import gld.sim.SimModel;

import java.awt.*;
import java.awt.event.*;
import java.util.*;

/**
*
* Extension of StatisticsView showing the statistics in bar charts.
*
* @author Group GUI
* @version 1.0
*/

public class StatsBarView extends StatisticsView
{
	public StatsBarView(StatisticsController parent, StatisticsModel stats)
	{
		super(parent, stats);
	}
	
	protected void paintStats(Graphics g)
	{
	}

	protected void paintAreaChanged() {}
}