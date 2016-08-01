
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

package gld.utils;

import java.awt.*;
import java.awt.geom.*;

/**
 * @author Pepijn van Lammeren
 * @version 1.0
 *
 * An extension of QuadCurve to draw roads, and calculate roaduser positions
 */
public class TurnCurve extends QuadCurve2D.Double
{
	/** The PathIterator for this TurnCurve */
	protected PathIterator iterator;
	/** The start point of this TurnCurve */
	protected Point start;
	/** The control point of this TurnCurve */
	protected Point control;
	/** The end point of this TurnCurve */
	protected Point end;
	/** A flag indicating whether this TurnCurve is actually a straight line */
	protected boolean straight;

	/**
	 * Create a new TurnCurve with specified end and control points
	 * @param p1 The start point
	 * @param cp The control point
	 * @param p2 The end point
	 */
	public TurnCurve(Point p1, Point cp, Point p2)
	{
		super((double)p1.x,(double)p1.y,(double)cp.x,(double)cp.y,(double)p2.x,(double)p2.y);
		start = p1;
		control = cp;
		end = p2;
		if(control.equals(avg(start,end)))
			straight = true;
		else
			straight = false;
		rewind();
	}
	
	/**
	 * Create a new straight TurnCurve with specified end points
	 * @param p1 The start point
	 * @param p2 The end point
	 */
	public TurnCurve(Point p1, Point p2)
	{
		this(p1,avg(p1,p2),p2);
	}
	
	/** Return the first point on this TurnCurve */
	public Point getFirst() { return start; }
	/** Return the control point of this TurnCurve */
	public Point getControl() { return control; }
	/** Return the last point on this TurnCurve */
	public Point getLast() { return end; }
	/** Returns true if this TurnCurve is actually a straight line, false otherwise */
	public boolean isStraight() { return straight; }

	/** Set the PathIterator to the first position on this TurnCurve
	 */
	public void rewind()
	{
		iterator = getPathIterator(null,0.001);
	}
	
	/** Let the PathIterator do a step and if successfull return the coordinates of it's current position
	 * @return The coordinats of the PathIterator's position, or null if no step was possible
	 */
	public Point next()
	{
		double[] c = new double[6];
		iterator.currentSegment(c);
		if(Math.round(c[0])==end.x && Math.round(c[1])==end.y)	// no more points
			return null;
		iterator.next();
		iterator.currentSegment(c);
		return new Point((int)Math.round(c[0]),(int)Math.round(c[1]));
	}

	
	/**
	 * Return the coordinates of the point with a given index on this TurnCurve
	 */
	public Point pointAt(int index)
	{
		rewind();
		for(int i=0; i < index; i++)
			iterator.next();
		double[] c = new double[6];
		int t = iterator.currentSegment(c);
		return new Point((int)Math.round(c[0]),(int)Math.round(c[1]));
	}


	private static Point avg(Point p, Point q)
	{
		return new Point((int)((p.x + q.x) / 2),(int)((p.y + q.y) / 2));
	}	
}