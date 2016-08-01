
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

import gld.sim.SimModel;
import gld.infra.Infrastructure;

import java.awt.Canvas;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.Color;
import java.awt.image.BufferedImage;
import java.awt.event.*;
import java.io.*;
import java.util.Observer;
import java.util.Observable;

/**
 *
 * A <code>TrackingView</code> shows a certain tracking graph.
 *
 * @author  Group GUI
 * @version 1.0
 */

public abstract class TrackingView extends Canvas implements Observer
{
	/** The separator used when saving data. */
	public static String SEP = "\t";
	/** Determines how many samples are tracked. */
	protected final static int MAX_DATA = 1000;
	/** Determines how many different data flows can be tracked. */
	protected int MAX_TRACK = getMaxTrack();
	/** Colors used to draw the graph. */
	protected Color[] colors = getColors();
	/** Determines whether a graph is drawn or not. */
	protected boolean[] show;

	/** The back buffer */
	protected BufferedImage buffer;
	/** The history of tracking data */
	protected float[][] trackData;
	/** The cycle this view started tracking */
	protected int startCycle;
	/** The cycle the end of the x-axis currently represents */
	protected int endCycle;
	  
	/** Used to determine where the user clicked the mouse. */
	protected float xStep;
	/** Index to the point in trackData to highlight. */
	protected int highlightIndex;
	/** Strings to draw in upper right corner. */
	protected String[] highlightText;
	/** String containing highlighted cycle text */
	protected String highlightCycle;
	  
	/** Current index in trackData[][]. */
	protected int curIndex;
	/** Highest valued sampled. */
	protected float maximum;
	/** Sample rate: each trackRate cycles, a sample is taken. */
	protected int trackRate;
	/** Current count: sample if trackCount == trackRate. */
	protected int trackCount;
  
	/**
	 * Creates a new <code>TrackingView</code>. A <code>TrackingView</code> is a component
	 * displaying a tracking graph.
	 *
	 * @param _startCycle The cycle this view starts tracking.
	 */
	public TrackingView(int _startCycle)
	{
		buffer = new BufferedImage(1,1,BufferedImage.TYPE_INT_BGR);

		trackData = new float[MAX_TRACK][MAX_DATA];
		highlightText = new String[MAX_TRACK];
		addMouseListener(new MyMouseListener());
		
		show = new boolean[MAX_TRACK];
		show[0] = true;
		for(int i=1; i<MAX_TRACK; i++)
			show[i] = false;
		
		init(_startCycle);
	}

	protected void init(int _startCycle)
	{
		startCycle = _startCycle;
		trackRate = 1;
		trackCount = 0;
		curIndex = 0;
		maximum = (float)5.0;
		endCycle = startCycle + MAX_DATA;
		highlightIndex = 0;
	}
	
	public void reset() { init(0); }




	/*============================================*/
	/* Basic GET and SET methods                  */
	/*============================================*/

	/** Returns the cycle this window started tracking */
	public int getStartCycle() { return startCycle; }
	/** Sets the cycle this window started tracking */
	public void setStartCycle(int cycle) { startCycle = cycle; }

	/** Determines the number of tracking graphs. */
	protected int getMaxTrack() { return 1; }
	/** Determines the colors of the tracking graphs. */
	protected Color[] getColors() { Color[] c = {Color.black}; return c; }
	/** Either shows or hides tracking graph <code>index</code>. */
	public void showGraph(int index, boolean b) { show[index] = b; }






	/*============================================*/
	/* Abstract methods                           */
	/*============================================*/

	/** Returns the next sample to be 'tracked'. */
	protected abstract float nextSample(int sourceIndex);

	/** Returns a short description of the specified source. */
	protected abstract String getSourceDesc(int sourceIndex);
	
	/** Returns the description for this tracking window. */
	public abstract String getDescription();
	
	/** Returns the Y-axis label */
	protected abstract String getYLabel();






	/*============================================*/
	/* Painting                                   */
	/*============================================*/


	/**
	 * This method is called by the <code>notifyObservers()</code> method
	 * in the <code>Model</code> representing the infrastructure.
	 *
	 * @param o The observed <code>Model</code>.
	 * @param arg Optional arguments.
	 */
	public void update(Observable o, Object arg) 
	{
		if(++trackCount != trackRate) return;
		trackCount = 0;
		
		for(int i=0; i<MAX_TRACK; i++) {
			trackData[i][curIndex] = nextSample(i);
			maximum = Math.max(maximum, trackData[i][curIndex]);
		}
		if(++curIndex == MAX_DATA)
			collapseData();
		
		redraw();
	}
	
	/** Redraws the buffer and paints it to the screen. */
	public void redraw() 
	{
		if (!isVisible()) return;
		
		drawBuffer();
		repaint();
	}

	/** Draws the buffer. */
	public void drawBuffer()
	{
		int w = buffer.getWidth();
		int h = buffer.getHeight();
		Graphics g = buffer.getGraphics();
		g.setColor(Color.white);
		g.fillRect(0,0,w,h);

		xStep = (float)(w - 40) / MAX_DATA;
		float yStep = (float)(h - 35) / maximum;

		g.setColor(Color.black);

		// y-axis
		String max = maximum + "";
		int dot;
		if((dot = max.indexOf(".")) > -1)
			max = max.substring(0, dot+2);
		g.drawString(max, 2, 15);
		g.drawString("0", 16, h - 22);
		g.drawString(getYLabel(), 0, 28);
		g.drawLine(30, 5, 30, h - 20);
		
		// x-axis
		g.drawString(startCycle + "", 30, h - 7);
		g.drawString(endCycle + "", w - 40, h - 7);
		g.drawString("time (cycles)", w / 2 - 20, h - 5);
		g.drawLine(30, h - 19, w - 10, h - 19);

		// the data
		for(int num=0; num<MAX_TRACK; num++)
			if(show[num]) {
				g.setColor(colors[num]);		
				for(int i=0; i<curIndex; i++)
					g.fillArc((int)(i * xStep) + 30, h - 20 - (int)(trackData[num][i] * yStep), 2, 2, 0, 360);
			}

		// highlight a certain point and give its values in trackData
		if(highlightIndex != 0)
		{
			g.setColor(Color.cyan);
			g.drawLine((int)(highlightIndex * xStep) + 30, 0, (int)(highlightIndex * xStep) + 30, h);
			g.setColor(Color.black);
			g.drawString(highlightCycle, w - 80, 15);
			int y = 35;
			for(int num=0; num<MAX_TRACK; num++)
				if(show[num]) {
					g.setColor(colors[num]);
					g.drawString(highlightText[num], w - 80, y);
					y += 15;
				}
		}
	}
	
	/** Override default update() to avoid flickering. */
	public void update(Graphics g)
	{
		paint(g);
	}

	/**
	 * Draws the data in the graphics buffer onto the screen.
	 *
	 * @param g the <code>Graphics</code> objects of this view.
	 */
	public void paint(Graphics g)
	{
		// draw buffer
		while (!g.drawImage(buffer,0,0,this));
	}
















	/*============================================*/
	/* Miscellaneous                              */
	/*============================================*/

	/**
	 * Sets the size for this <code>View</code>
	 *
	 * @param width the new width
	 * @param height the new height
	 */
	public void setSize(int width, int height)
	{
		buffer = new BufferedImage(width+1,height+1,BufferedImage.TYPE_INT_BGR);
		super.setSize(width,height);
		redraw();
	}

	public void setBounds(int x, int y, int width, int height)
	{
		width = width < 1 ? 1 : width;
		height = height < 1 ? 1 : height;
		buffer = new BufferedImage(width+1,height+1,BufferedImage.TYPE_INT_BGR);
		super.setBounds(x,y,width,height);
		redraw();
	}
	
	/**
	* Collapses trackData and changes other variables accordingly.
	*/
	protected void collapseData()
	{
		for(int i=0; i<MAX_TRACK; i++)
			for(int j=0; j<MAX_DATA/2; j++)
				trackData[i][j] = trackData[i][j*2];
		curIndex = MAX_DATA / 2;
		trackRate *= 2;
		endCycle = startCycle + (endCycle - startCycle) * 2;
		highlightIndex = 0;
	}
	
	/**
	 * Save data to a CSV file. 
	 */
	public void saveData(String filename, SimModel model) throws IOException
	{
		PrintWriter out=new PrintWriter(new FileWriter(new File(filename)));
		out.println("# Data exported by Green Light District"); out.println("#");
		Infrastructure infra = model.getInfrastructure();
		out.println("# Infrastructure: \"" + infra.getTitle() + "\" by " + infra.getAuthor());
		out.println("# Simulation: \"" + model.getSimName() + "\"");
		out.println("# Tracking data: \"" + getDescription() + "\""); 
		out.println("# Measured: " + getYLabel()); out.println("#");
		String fieldsDesc = "";
		for(int i=0; i<MAX_TRACK; i++)
			if(show[i])
				fieldsDesc += "<" + getSourceDesc(i) + ">" + (i==MAX_TRACK-1?"":SEP);
		out.println("# Format: <Cycle>" + SEP + fieldsDesc); out.println("#");
		
		String ln;
		for(int i=0; i<curIndex; i++) {
			ln = (startCycle + i * trackRate) + SEP;
			for(int j=0; j<MAX_TRACK; j++)
				if(show[j])
					ln += trackData[j][i] + (j==MAX_TRACK-1?"":SEP);
			out.println(ln);
		}
		out.close();
	}






	/*============================================*/
	/* Listeners                                  */
	/*============================================*/
	
	/** 
	* Listens to mouse clicks on this TrackingView. 
	* Sets the index of the graph datatable to be highlighted.
	*/
	protected class MyMouseListener extends MouseAdapter
	{
		public void mouseClicked(MouseEvent e)
		{
			if((e.getModifiers() & InputEvent.BUTTON3_MASK) > 0) // right-click
			{
				highlightIndex = 0;
			}
			else
			{
				Point p = e.getPoint();
				highlightIndex = (int)((p.x - 30) / xStep);
				if(highlightIndex < 0 || highlightIndex > curIndex)
				{
					highlightIndex = 0;
					return;
				}
				highlightCycle = "cycle = " + (startCycle + highlightIndex * trackRate);
				String desc; 
				for(int i=0; i<MAX_TRACK; i++) {
					highlightText[i] = getSourceDesc(i) + " = " + trackData[i][highlightIndex];
					int dot;
					if((dot = highlightText[i].indexOf(".")) > -1)
						highlightText[i] = highlightText[i].substring(0,dot+2);
				}
			}
			redraw();
		}
	}
}