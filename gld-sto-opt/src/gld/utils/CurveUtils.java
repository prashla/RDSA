
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
import java.util.*;
import gld.infra.InfraException;
import gld.infra.Road;
import gld.infra.Node;
import gld.infra.Turn;
import gld.utils.CurveException;
import gld.utils.TurnCurve;

/**
 * This class presents methods to handle curves
 * @author Pepijn van Lammeren
 * @version 1.0
 */
public class CurveUtils
{
	/**
	 * Calculate the intersection coordinats for two lines, given a base point and angle for
	 * each line.
	 * @param p1 The base point for the first line
	 * @param p2 The base point for the second line
	 * @param ar1 The angle of the first line in radians
	 * @param ar2 The angle of the second line in radians
	 */
	public static Point calcIntersect(Point p1, Point p2, double ar1, double ar2) throws CurveException
	{
   		Point p12 = new Point((int)(p1.x + 1000 * Math.cos(ar1)),(int)(p1.y - 1000 * Math.sin(ar1)));
   		Point p22 = new Point((int)(p2.x + 1000 * Math.cos(ar2)),(int)(p2.y - 1000 * Math.sin(ar2)));
   		double m1,m2,b1,b2,px,py;
   		if(p1.x == p12.x && p2.x == p22.x)
			return new Point((int)((p1.x + p2.x) / 2) + 2,(int)((p1.y + p2.y) / 2) + 2);
   		else if(p1.x == p12.x)
	   	{
			m2 = (double)(p22.y - p2.y) / (p22.x - p2.x);
			b2 = (double)p2.y - (m2 * p2.x);
			px = (double) p12.x;
			py = (double) m2 * px + b2;
		}
		else if(p2.x == p22.x)
		{
			m1 = (double)(p12.y - p1.y) / (p12.x - p1.x);
			b1 = (double) p1.y - (m1 * p1.x);
			px = (double) p12.x;
			py = (double) m1 * px + b1;
		}
		else
		{
			m1 = (double)(p12.y - p1.y) / (p12.x - p1.x);
			m2 = (double)(p22.y - p2.y) / (p22.x - p2.x);
			if(m1==m2)
				return new Point((int)((p1.x + p2.x) / 2),(int)((p1.y + p2.y) / 2));
			b1 = (double)p1.y - (m1 * p1.x);
			b2 = (double)p2.y - (m2 * p2.x);
			px = (b2 - b1) / (m1 - m2);
			py = m1 * px + b1;
		}
		return new Point((int)px,(int)py);
	}


	public static TurnCurve createCurve(Point p1, Point p2, double a1, double a2) throws CurveException
	{
		Point tp = new Point((int)(p1.x + p1.distance(p2) * Math.cos(a1)),(int)(p1.y - p1.distance(p2) * Math.sin(a1)));
		Rectangle r = new Rectangle(p2.x - 10, p2.y - 10, 20, 20);
		if(r.contains(tp))
			return new TurnCurve(p1,p2);
		Point ip = calcIntersect(p1,p2,a1,a2);
		return new TurnCurve(p1,ip,p2);
	}

	/**
	 * Creates a QCurve given two points and angles
	 * @param p1 The coordinats of the first point
	 * @param p2 The coordinats of the second point
	 * @param a1 The angle of the first point in radians
	 * @param a2 The angle of the second point
	 *
	public static QCurve createCurve(Point p1, Point p2, double a1, double a2) throws CurveException
	{
		Point tp = new Point((int)(p1.x + p1.distance(p2) * Math.cos(a1)),(int)(p1.y - p1.distance(p2) * Math.sin(a1)));
		Rectangle r = new Rectangle(p2.x - 10, p2.y - 10, 20, 20);
		Point ip;
		if(r.contains(tp))
			return new QCurve(p1,p2);
		ip = calcIntersect(p1,p2,a1,a2);
		QCurve c = new QCurve(p1,ip,p2);
// test voor path-iterator-bug:
		Point2D.Double p = c.next(c.getFirst(),2);
		if(!(c.getBounds().contains(p)))
		{
			c = new QCurve(p2,ip,p1);
			p = c.next(c.getFirst(),2);
			c.setSwitched(true);
		}
		return c;
	}

	
	/**
	 * Creates a GeneralPath connecting both curves in a Drivelane-like fashion
	 * @param c1 The first curve
	 * @param c2 The second curve
	 * @return A GeneralPath connecting both curves
   */
	public static GeneralPath createPath(TurnCurve c1,TurnCurve c2)
	{
		GeneralPath p = new GeneralPath(c1);
		p.lineTo((float)c2.getX2(),(float)c2.getY2());
		p.quadTo((float)c2.getCtrlX(),(float)c2.getCtrlY(),(float)c2.getX1(),(float)c2.getY1());
		p.closePath();
		return p;
	}
		
	
	
	/**
	 * Calculates the length of the specified curve in steps of size step
	 * @param c The QCurve whose length to calculate
	 * @param step The stepsize
	 * @return The length of the specified QCurve
	 */
	public static double calcLength(TurnCurve c, int step)
	{
		Vector points = new Vector();
		c.rewind();
		Point p;
		points.add(c.getFirst());
		while((p = c.next()) != null)
			points.add(p);
		double d = 0.0;
		for(int i = 1; i < points.size(); i++)
		{
			Point p1 = (Point)(points.get(i-1));
			Point p2 = (Point)(points.get(i));
			d += p2.distance(p1);
		}
		return (d/(double)step);
	}
	
	
	/**
	 * Returns the number of steps used to traverse the specified curve
	 * @param c The QCurve
	 * @return The number of steps
	 */
	public static int calcPathSteps(TurnCurve c)
	{
		if(c.isStraight())
			return (int)(c.getFirst().distance(c.getLast()));
		c.rewind();
		int n = 0;
		Point p;
		while((p = c.next()) != null)
			n++;
		return n;
	}
	
    
   /**
    * Returns the coordinates of the point with specified index and stepsize
    * on a given curve
    * @param c The curve
    * @param step The stepsize
    * @param index The index of the point
    * @return The Point at the specified index
   */
	public static Point getPoint(TurnCurve c, int index, int step)
	{
		if(c.isStraight())
		{
			double dx = (double)(c.getLast().x - c.getFirst().x) / (double)((double)c.getFirst().distance(c.getLast()) / (double)step);
			double dy = (double)(c.getLast().y - c.getFirst().y) / (double)((double)c.getFirst().distance(c.getLast()) / (double)step);
			Point p = new Point((int)((double)c.getFirst().x + dx * (double)index),
				  	            (int)((double)c.getFirst().y + dy * (double)index));
			if(c.getFirst().distance(p) > c.getFirst().distance(c.getLast()))
				return c.getLast();
			if(c.getLast().distance(p) > c.getLast().distance(c.getFirst()))
				return c.getFirst();			
			return p;
		}
		Point p;
		Point q = c.getFirst();
		double pos = 0.0;
		while((p = c.next()) != null)
		{
			pos += p.distance(q);
			if(pos >= index * step)
				return p;
			q = p;
		}
		return c.getLast();
 	}

    /**
     * Returns the angle of a given curve at the given index
     * @param c The curve
     * @param step The stepsize
     * @param index The index
     * @return The angle of the curve at the specified index
    */
    public static double getAngle(TurnCurve c, int index, int step)
    {
		Point p1;
		Point p2;
		if(c.isStraight())
		{
			p1 = c.getFirst();
			p2 = c.getLast();
		}
		else
		{			
			p1 = getPoint(c,index,step);
			p2 = null;
			for(int n = 0; n < 10; n++)
				p2 = c.next();
		
			if(p2==null)
			{
				p2 = p1;
				p1 = c.getControl();
			}			
		}
		if(p2.y - p1.y == 0)
		{
			if(p1.x > p2.x)
				return (double)(-Math.PI * 0.5);
			else
				return (double)(Math.PI * 0.5);
		}
		double a = Math.atan((double)(p2.x - p1.x) / (double)(p2.y - p1.y));
		return a;
	}
		


	public static double getCurrentAngle(TurnCurve c, Point p)
	{
		Point p1;
		Point p2;
		if(c.isStraight())
		{
			p1 = c.getFirst();
			p2 = c.getLast();
		}
		else
		{
			p1 = p;
			p2 = null;
			for(int n = 0; n < 10; n++)
				p2 = c.next();
			
			if(p2==null)
			{
				p2 = p1;
				p1 = c.getControl();
			}
		}
		if(p2.y - p1.y == 0)
		{
			if(p1.x > p2.x)
				return (double)(-Math.PI * 0.5);
			else
				return (double)(Math.PI * 0.5);
		}
		return Math.atan((double)(p2.x - p1.x) / (double)(p2.y - p1.y));
	}


    /**
     * Convert any angle to a value in the range [0, 2Pi]
     */
	public static double normalize(double a)
	{
		return a > 0 ? a % (Math.PI * 2) : (a % (Math.PI * 2)) + Math.PI * 2;
	}



	public static void setupRoadSizes(Road road, int step) throws InfraException
	{
		int acp, bcp;
    	acp = road.getAlphaNode().isConnectedAt(road);
    	bcp = road.getBetaNode().isConnectedAt(road);

		Node alphaNode = road.getAlphaNode();
		Node betaNode = road.getBetaNode();

		int alphaMaxWidth = alphaNode.getWidth();
		int betaMaxWidth = betaNode.getWidth();

		int cax=0,cay=0,cbx=0,cby=0;
		double alphaAngle, betaAngle;

	// calculate the starting positions of the road
		if(acp == 0)    // alpha top
		{
        	cax = alphaNode.getCoord().x;
        	cay = alphaNode.getCoord().y - alphaMaxWidth * 5 - 1;
        	alphaAngle = Math.toRadians(90);
    	}
    	else if(acp == 1) // alpha right
    	{
            cax = alphaNode.getCoord().x + alphaMaxWidth * 5 + 1;
            cay = alphaNode.getCoord().y;
            alphaAngle = 0.0;
        }
        else if(acp == 2) // alpha bottom
        {
            cax = alphaNode.getCoord().x;
            cay = alphaNode.getCoord().y + alphaMaxWidth * 5 + 1;
            alphaAngle = Math.toRadians(270);
        }
        else    // alpha left
        {
            cax = alphaNode.getCoord().x - alphaMaxWidth * 5 - 1;
            cay = alphaNode.getCoord().y;
            alphaAngle = Math.toRadians(180);
        }
		
        if(bcp == 0)    // beta top
        {
            cbx = betaNode.getCoord().x;
            cby = betaNode.getCoord().y - betaMaxWidth * 5 - 1;
            betaAngle = Math.toRadians(270);
        }
        else if(bcp == 1) // beta right
        {
            cbx = betaNode.getCoord().x + betaMaxWidth * 5 + 1;
            cby = betaNode.getCoord().y;
            betaAngle = Math.toRadians(180);
        }
        else if(bcp == 2) // beta bottom
        {
            cbx = betaNode.getCoord().x;
            cby = betaNode.getCoord().y + betaMaxWidth * 5 + 1;
            betaAngle = Math.toRadians(90);
        }
        else    // beta left
        {
            cbx = betaNode.getCoord().x - betaMaxWidth * 5 - 1;
            cby = betaNode.getCoord().y;
            betaAngle = 0.0;
        }

        Point alphaCoord = new Point(cax,cay);
        Point betaCoord = new Point(cbx,cby);
	
		if(road.getTurns().length==0)	// no turns
		{
			road.setLength((int)(alphaCoord.distance(betaCoord) / step));
			return;
		}
		Turn[] turns = road.getTurns();
		double pos = 0.0;

		try
		{
			TurnCurve arc = createCurve(alphaCoord,turns[0].getCoord(),alphaAngle,turns[0].getAngle() + Math.PI);
			pos = calcLength(arc,step);
			turns[0].setPosition((int)pos);
		}
		catch(CurveException ce)
		{
		}
		
		for(int i = 0; i < turns.length - 1; i++)
		{
			try
			{
				TurnCurve arc = createCurve(turns[i].getCoord(),turns[i+1].getCoord(),turns[i].getAngle(),turns[i+1].getAngle() + Math.PI);
				pos += calcLength(arc,step);
				turns[i+1].setPosition((int)pos);
			}
			catch(CurveException ce)
			{
			}
		}
	
		try
		{
			TurnCurve arc = createCurve(turns[turns.length-1].getCoord(),betaCoord,turns[turns.length-1].getAngle(),betaAngle);
			pos += calcLength(arc,step);
		}
		catch(CurveException ce)
		{
		}

		road.setLength((int)Math.round(pos));
	}			

    /**
     * Returns the length of a given Road.
     * @param road A Road
     * @param step The step-size in pixels
     */
	    public static int getRoadLength(Road road, int step) throws InfraException
    	{	

		if(road.getTurns().length==0)	// no turns
		{
	    	int dx = (road.getAlphaNode().getCoord().x - road.getBetaNode().getCoord().x);
		    int dy = (road.getAlphaNode().getCoord().y - road.getBetaNode().getCoord().y);
		    return (int)(Math.sqrt(dx*dx+dy*dy) / step);
		}
		Turn[] turns = road.getTurns();
		double len = 0;

		int acp, bcp;
    	acp = road.getAlphaNode().isConnectedAt(road);
    	bcp = road.getBetaNode().isConnectedAt(road);

		Node alphaNode = road.getAlphaNode();
		Node betaNode = road.getBetaNode();

		int alphaMaxWidth = alphaNode.getWidth();
		int betaMaxWidth = betaNode.getWidth();

		int cax=0,cay=0,cbx=0,cby=0;
		double alphaAngle, betaAngle;

	// calculate the starting positions of the road
		if(acp == 0)    // alpha top
		{
        	cax = alphaNode.getCoord().x;
        	cay = alphaNode.getCoord().y - alphaMaxWidth * 5 - 1;
        	alphaAngle = Math.toRadians(90);
    	}
    	else if(acp == 1) // alpha right
    	{
            cax = alphaNode.getCoord().x + alphaMaxWidth * 5 + 1;
            cay = alphaNode.getCoord().y;
            alphaAngle = 0.0;
        }
        else if(acp == 2) // alpha bottom
        {
            cax = alphaNode.getCoord().x;
            cay = alphaNode.getCoord().y + alphaMaxWidth * 5 + 1;
            alphaAngle = Math.toRadians(270);
        }
        else    // alpha left
        {
            cax = alphaNode.getCoord().x - alphaMaxWidth * 5 - 1;
            cay = alphaNode.getCoord().y;
            alphaAngle = Math.toRadians(180);
        }
		
        if(bcp == 0)    // beta top
        {
            cbx = betaNode.getCoord().x;
            cby = betaNode.getCoord().y - betaMaxWidth * 5 - 1;
            betaAngle = Math.toRadians(270);
        }
        else if(bcp == 1) // beta right
        {
            cbx = betaNode.getCoord().x + betaMaxWidth * 5 + 1;
            cby = betaNode.getCoord().y;
            betaAngle = Math.toRadians(180);
        }
        else if(bcp == 2) // beta bottom
        {
            cbx = betaNode.getCoord().x;
            cby = betaNode.getCoord().y + betaMaxWidth * 5 + 1;
            betaAngle = Math.toRadians(90);
        }
        else    // beta left
        {
            cbx = betaNode.getCoord().x - betaMaxWidth * 5 - 1;
            cby = betaNode.getCoord().y;
            betaAngle = 0.0;
        }

        Point alphaCoord = new Point(cax,cay);
        Point betaCoord = new Point(cbx,cby);

	
       
        try
        {
            TurnCurve arc = CurveUtils.createCurve(alphaCoord,turns[0].getCoord(),alphaAngle,normalize(turns[0].getAngle()));
		    len = len + calcLength(arc,step);
        }
        catch(CurveException exc)
        {
			System.out.println("CurveException");
	    	int dx = alphaCoord.x - turns[0].getCoord().x;
	    	int dy = alphaCoord.y - turns[0].getCoord().y;
	    	len = len + (Math.sqrt(dx*dx+dy*dy) / step);
        }
		for(int t = 0; t < turns.length - 1; t++)
    	{
            double ta1 = CurveUtils.normalize(turns[t].getAngle());
            double ta2 = CurveUtils.normalize(turns[t+1].getAngle() + Math.PI);
            Point tp1 = turns[t].getCoord();
            Point tp2 = turns[t+1].getCoord();
            try
            {
                TurnCurve arc = createCurve(tp1,tp2,ta1,ta2);
				len = len + calcLength(arc,step);
            }
            catch(CurveException exc)
            {
				System.out.println("CurveException");
				int dx = tp1.x - tp2.x;
				int dy = tp1.y - tp2.y;
				len = len + (Math.sqrt(dx*dx+dy*dy) / step);
            }
		}
        try
        {
            TurnCurve arc = CurveUtils.createCurve(turns[turns.length-1].getCoord(),betaCoord,turns[turns.length-1].getAngle(),betaAngle);
  	    	len = len + calcLength(arc,step);
        }
        catch(CurveException exc)
        {
			System.out.println("CurveException");	    	
	    	int dx = betaCoord.x - turns[turns.length-1].getCoord().x;
	    	int dy = betaCoord.y - turns[turns.length-1].getCoord().y;
	    	len = len + (int)(Math.sqrt(dx*dx+dy*dy)/step);
        }

		return (int)len;
	}
}