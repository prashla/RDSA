
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
* The view of the statistics viewer, controlled by StatisticsController.
*
* @author Group GUI
* @version 1.0
*/

public abstract class StatisticsView extends Canvas implements Observer
{
	/** Separator used when saving data. */
	public static String SEP = "\t";
	
	protected StatisticsController parent;
	protected StatisticsModel stats;

	protected int horScroll, verScroll;
	protected Rectangle paintArea;

	protected final Font titleFont = new Font("arial",Font.BOLD,14);
	protected final Font infoFont = new Font("arial",Font.PLAIN,12);
	protected final Font tableFont = new Font("arial",Font.PLAIN,11);


	/**
	* Creates a <code>StatisticsView</code>.
	*
	* @param _model The <code>SimModel</code> statistics should be read from.
	*/
	public StatisticsView(StatisticsController _parent, StatisticsModel _stats)
	{
		parent = _parent;
		stats = _stats;
		horScroll = verScroll = 0;
		paintArea = getBounds();
	}


	/** Invoked when the StatisticsModel is changed. */
	public void update(Observable obs, Object obj) { update(); }
	
	/** Invoked when the view should be redone. */
	public void update()
	{
		paintAreaChanged();
		repaint();
	}
	
	/** Paints the view. */
	public void paint(Graphics g) 
	{
		g.setColor(Color.white);
		g.fillRect(0,0,getWidth(),getHeight());
		g.setColor(Color.black);
		g.setFont(titleFont);
		g.drawString("Statistics for simulation \"" + stats.getSimName() + "\" (at cycle " + stats.getCycle() + ")", 20, 20);
		g.setFont(infoFont);
		g.drawString("Infrastructure: \"" + stats.getInfraName() + "\" by " + stats.getInfraAuthor(), 20, 40);
	
		paintStats(g);
	}




	/*============================================*/
	/* GET AND SET                                */
	/*============================================*/

	/** Sets the StatisticsModel to be shown. */
	public void setStatisticsModel(StatisticsModel _stats) { stats = _stats; update(); }
	/** Returns the StatisticsModel to be shown. */
	public StatisticsModel getStatisticsModel() { return stats; }
	/** Sets horScroll, the horizontal scrolling value. */
	public void setHorScroll(int hs) { horScroll = hs; update(); }
	/** Returns horScroll, the horizontal scrolling value. */
	public int getHorScroll() { return horScroll; }
	/** Sets verScroll, the vertical scrolling value. */
	public void setVerScroll(int vs) { verScroll = vs; update(); }
	/** Returns verScroll, the vertical scrolling value. */
	public int getVerScroll() { return verScroll; }
	/** Sets the extension-specific paint area of this view. */
	protected void setPaintArea(Rectangle r) { paintArea = r; update(); }
	/** Returns the extension-specific paint area of this view. */
	protected Rectangle getPaintArea() { return paintArea; }
	/** Overrides default method to change paintArea accordingly. */
	public void setBounds(int x, int y, int w, int h)
	{
		int paw = Math.max(1, w - 40), pah = Math.max(1, h - 75);
		setPaintArea(new Rectangle(20, 55, paw, pah));
		super.setBounds(x,y,w,h);
	}





	/*============================================*/
	/* ABSTRACT METHODS                           */
	/*============================================*/

	/** Invoked when the paintArea is changed. */
	protected abstract void paintAreaChanged();

	/** Paints the statistics in Rectangle r on the view. */
	protected abstract void paintStats(Graphics g);
}