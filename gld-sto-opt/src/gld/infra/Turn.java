
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

import gld.xml.*;
import java.awt.Point;
import java.io.IOException;
import java.util.*;

/**
 *
 * Basic turn.
 * A turn may later be sophisticated with bends in the road.
 *
 * @author Group Datastructures
 * @version 1.0
 */

public class Turn implements XMLSerializable
{
	/** The coordinates of this turn */
	protected Point point;
	/** The angle of this turn */
	protected double angle;
	protected String parentName="model.infrastructure.node.road";

	protected int type;
	/** The position (index) of this turn on the road */
	protected double position;

	public static final int TURNPOINT = 1;
	public static final int CONTROLPOINT = 2;

	
	public Turn(Point p, double a) {
		point = p;
		angle = a;
	}
	
	public Turn(Point p, double a, int t) {
		point = p;
		angle = a;
		type = 1;
	}
	
	public Turn(Point p, int t)
	{
		point = p;
		type = t;
	}
	
	public Turn() { }

	/** Returns the coordinates of this turn */
	public Point getCoord() { return point; }
	/** Sets the coordinates of this turn */
	public void setCoord(Point p) { point = p; }

	/** Returns the angle of this turn */
	public double getAngle() { return angle; }
	/** Sets the angle of this turn */
	public void setAngle(double a) { angle = a; }

	/** Returns the type of this turn */
	public int getType() { return type; }

	/** Returns the position of this turn on the road */
	public double getPosition() { return position; }
	/** Sets the position of this turn on the road */
	public void setPosition(double p) { position = p; }

  // XML Serializable implementation
  
	public void load (XMLElement myElement,XMLLoader loader) throws XMLTreeException,IOException,XMLInvalidInputException
	{	point=new Point (myElement.getAttribute("x-pos").getIntValue(),
                		 myElement.getAttribute("y-pos").getIntValue());
		angle=myElement.getAttribute("angle").getDoubleValue();		 
		position=myElement.getAttribute("rel-pos").getDoubleValue();
	}

	public XMLElement saveSelf () throws XMLCannotSaveException
	{	XMLElement result=new XMLElement("turn");
		result.addAttribute(new XMLAttribute("x-pos", (int)(point.getX())));
		result.addAttribute(new XMLAttribute("y-pos", (int)(point.getY())));
		result.addAttribute(new XMLAttribute("angle",angle));
		result.addAttribute(new XMLAttribute("rel-pos", position));
		return result;
	}
 
	public void saveChilds (XMLSaver saver) throws XMLTreeException,IOException,XMLCannotSaveException
	{ 	// A turn has no childs
	}

	public String getXMLName () 
	{	return parentName+".turn";
	}
	
	public void setParentName (String parentName)
	{	this.parentName=parentName; 
	}
	
}
