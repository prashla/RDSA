
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

package gld.infra;

import gld.*;
import gld.xml.*;
import java.awt.Graphics;
import java.awt.Color;
import java.io.IOException;
import java.util.*;

/**
 * A Car is the standard menace on our Roads.
 * 
 * @author Group Datastructures
 * @version 1.0
 *
 * A Car is the standard menace on our Roads.
 */

public class Car extends Automobile
{
	protected final int type = RoaduserFactory.getTypeByDesc("Car");
	protected final int length = 2;
	protected final int passengers = 1;
	protected int speed = 2;


	public Car(Node new_startNode, Node new_destNode, int pos) {
		super(new_startNode, new_destNode, pos);
		color = RoaduserFactory.getColorByType(type);
		int r = color.getRed();
		int g = color.getGreen();
		int b = color.getBlue();
		if(r==0) r = (int)(Math.random() * 160);
		if(g==0) g = (int)(Math.random() * 160);
		if(b==0) b = (int)(Math.random() * 160);
		color = new Color(r,g,b);
	}
	
	/** Empty constructor for loading */
	public Car() { }
	
	public String getName() { return "Car"; }
	
	/** Returns the speed of this Roaduser in blocks per cycle */
	public int getSpeed() { return speed; }
	public int getLength() { return length; }
	public int getType() { return type; }
	public int getNumPassengers() { return passengers; }	
	
	public void paint(Graphics g, int x, int y, float zf) {
		paint(g,x,y,zf,(double)0.0);
	}

	public void paint(Graphics g, int x, int y, float zf, double dlangle) {
		g.setColor(color);
		double angle = dlangle - Math.toRadians(45.0);
		int[] cx = new int[4];
		cx[0] = (int)(Math.round((double)x + Math.sin(angle) * 3));
		cx[1] = (int)(Math.round((double)x + Math.sin(angle + Math.toRadians(90.0)) * 3));
		cx[2] = (int)(Math.round((double)x + Math.sin(angle + Math.toRadians(180.0) + Math.atan(0.5)) * 6));
		cx[3] = (int)(Math.round((double)x + Math.sin(angle + Math.toRadians(270.0) - Math.atan(0.5)) * 6));
		
		int[] cy = new int[4];
		cy[0] = (int)(Math.round((double)y + Math.cos(angle) * 3));
		cy[1] = (int)(Math.round((double)y + Math.cos(angle + Math.toRadians(90.0)) * 3));
		cy[2] = (int)(Math.round((double)y + Math.cos(angle + Math.toRadians(180.0) + Math.atan(0.5)) * 6));
		cy[3] = (int)(Math.round((double)y + Math.cos(angle + Math.toRadians(270.0) - Math.atan(0.5)) * 6));
		
		g.fillPolygon(cx,cy,4);
		//g.fillRect((int)((x - 3) * zf),(int)((y - 3) * zf),(int) (7 * zf),(int) (7 * zf));
	}

	
    // Specific XMLSerializable implementation 
    
 	public String getXMLName ()
 	{ 	return parentName+".roaduser-car";
 	}
 
}