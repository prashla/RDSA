
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
 * The Bus. Red, Big, Mean and Lean people moving machine.
 *
 * @author Group Datastructures
 * @version 1.0
 *
 */

public class Bus extends Automobile
{
	protected final int type = RoaduserFactory.getTypeByDesc("Bus");
	protected final int length = 3;
	protected final int speed = 2;
	protected final int passengers = 10;

	public Bus(Node new_startNode, Node new_destNode, int pos) {
		super(new_startNode, new_destNode, pos);
		// make color little bit more random
		color = RoaduserFactory.getColorByType(type);
		int r = color.getRed();
		int g = color.getGreen();
		int b = color.getBlue();
		if(r==0) r = (int)(Math.random() * 160);
		if(g==0) g = (int)(Math.random() * 160);
		if(b==0) b = (int)(Math.random() * 160);
		color = new Color(r,g,b);
	}
	
	public Bus() { }
	
	public String getName() { return "Bus"; }
	
	/** Returns the speed of this Roaduser in blocks per cycle */
	public int getSpeed() { return speed; }
	public int getLength() { return length; }
	public int getType() { return type; }
	public int getNumPassengers() { return passengers; }	

	public void paint(Graphics g, int x, int y, float zf) {}
	
	public void paint(Graphics g, int x, int y, float zf, double dlangle) {
		g.setColor(color);
		double angle = dlangle - Math.toRadians(45.0);
		int[] cx = new int[4];
		cx[0] = (int)(Math.round((double)x + Math.sin(angle) * 3));
		cx[1] = (int)(Math.round((double)x + Math.sin(angle + Math.toRadians(90.0)) * 3));
		cx[2] = (int)(Math.round((double)x + Math.sin(angle + Math.toRadians(180.0) + Math.atan(1 - 1.0 / 3.0)) * 10));
		cx[3] = (int)(Math.round((double)x + Math.sin(angle + Math.toRadians(270.0) - Math.atan(1 - 1.0 / 3.0)) * 10));
		
		int[] cy = new int[4];
		cy[0] = (int)(Math.round((double)y + Math.cos(angle) * 3));
		cy[1] = (int)(Math.round((double)y + Math.cos(angle + Math.toRadians(90.0)) * 3));
		cy[2] = (int)(Math.round((double)y + Math.cos(angle + Math.toRadians(180.0) + Math.atan(1 - 1.0 / 3.0)) * 10));
		cy[3] = (int)(Math.round((double)y + Math.cos(angle + Math.toRadians(270.0) - Math.atan(1 - 1.0 / 3.0)) * 10));
		
		g.fillPolygon(cx,cy,4);
	}
	
	// Specific XMLSerializable implementation 
    
	public String getXMLName ()
	{ 	return parentName+".roaduser-bus";
	}
}