
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

import java.io.IOException;
import java.util.*;
import java.awt.Graphics;
import java.awt.*;
import java.applet.*;

/**
 *
 * Basic Roaduser
 *
 * @author Group Datastructures
 * @version 1.0
 *
 * Todo:
 * Fix Save/Load
 */

abstract public class Roaduser implements Selectable, XMLSerializable, TwoStageLoader, Cloneable
{
	/** The node this Roaduser spawned at */
	protected Node startNode;
	/** The node that is the destination of this Roaduser */
	protected Node destNode;
	/** The last cycle this Roaduser moved */
	protected int cycleMoved;
	/** The last cycle this Roaduser was asked when it had last moved */
	protected int cycleAsked;
	/** The position of this Roaduser on the drivelane. Zero based. */
	protected int position;
	/** The SignID of the lane this Roaduser came from */
	protected int prevSign = -1;
	/** The starttime on this lane */
	protected int drivelaneStartTime;
	/** The delay experienced so far */
	protected int delay;
	/** Stuff to transfer between the first and second stage loader */
	protected TwoStageLoaderData loadData=new TwoStageLoaderData();
	/** The color of this Roaduser */
	protected Color color = new Color(0, 0, 0);
	/** The name of the parent of this Roaduser */
	protected String parentName="model.infrastructure.lane";
	
    /** The id of the sign this roaduser last passed */
	protected int waitTl=-1;
	/** The id of the sign this roaduser last passed */
	protected int waitPos=-1;
	/** The id of the sign this roaduser last passed */
	protected boolean waitTlColor=false;
    /** The id of the sign this roaduser last passed */
	protected int prevWaitTl=-1;
	/** The id of the sign this roaduser last passed */
	protected int prevWaitPos=-1;
	/** The id of the sign this roaduser last passed */
	protected boolean prevWaitTlColor=false;
	protected boolean voted=false;
	protected boolean inQueueForSign=false;
	
	public void setInQueueForSign(boolean b) { inQueueForSign=b; }
	public boolean getInQueueForSign() { return inQueueForSign; }

	public Roaduser(Node _startNode, Node _destNode, int pos)
	{	this();
		startNode = _startNode;
		destNode = _destNode;
		position = pos;
	}
	
	public Roaduser() 
	{	resetStats(); 
	}
	
	public void resetStats ()
	{	cycleMoved = -1;
		cycleAsked = -1;
		drivelaneStartTime = -1;
		delay=0;
	}
	
	public Object clone () 
	{ 	try {return super.clone();} catch(Exception e) {System.out.println(e);}
		return null;
	}
		











	/*============================================*/
	/* Basic GET and SET methods                  */
	/*============================================*/

	/** Returns the start Node of this Roaduser */
	public Node getStartNode() { return startNode; }
	/** Sets the start Node of this Roaduser */
	public void setStartNode(Node n) { startNode = n; }

	/** Returns the destination Node of this Roaduser */
	public Node getDestNode() { return destNode; }
	/** Sets the destination Node of this Roaduser */
	public void setDestNode(Node n) { destNode = n; }

	/** Returns the position of this Roaduser on the current Drivelane */
	public int getPosition() { return position; }
	/** Sets the position of this Roaduser on the current Drivelane */
	public void setPosition(int pos) { /*System.out.println("newPos:"+pos);*/ position = pos; }
	
	/** Returns the last cycle this Roaduser moved */
	public int getCycleMoved() { return cycleMoved; }
	/** Sets the last cycle this Roaduser moved */
	public void setCycleMoved(int cycle) { cycleMoved = cycle; }
	
	/** Returns the last cycle this Roaduser was asked about its movements */
	public int getCycleAsked() { return cycleAsked; }
	/** Sets the last cycle this Roaduser was asked its movements */
	public void setCycleAsked(int cycle) { cycleAsked = cycle; }

	/** Returns the start time of this Roaduser in the current drivelane */
	public int getDrivelaneStartTime() { return drivelaneStartTime; }
	/** Sets the start time of this Roaduser on the current Drivelane */
	public void setDrivelaneStartTime(int time) { drivelaneStartTime = time; }

	/** Returns the distance experienced so far */
	public int getDelay() { return delay; }
	/** Add a given delay to the total delay already experienced */
	public void addDelay(int d) { delay += d; }
	/** Sets a new delay */
	public void setDelay(int delay) { this.delay=delay;}
	/** Returns the Id of the previous lane this Roaduser hit */
	public int getPrevSign() { return prevSign; }
	/** Sets the Id of the previous lane this Roaduser hit */
	public void setPrevSign(int _prevSign) { prevSign = _prevSign; }
	/** Returns the color of this Roaduser */
	public Color getColor(){ return color; }	
	/** Sets the color of this Roaduser */
	public void setColor(Color c) { color = c; }

	/** Sets the last waiting point's position relative to the sign*/
	public void setWaitPos(int tlId, boolean b, int pos) { 
	    prevWaitTl = waitTl;
	    prevWaitTlColor = waitTlColor;
	    prevWaitPos = waitPos;
	    waitTl=tlId; 
	    waitTlColor=b;
	    waitPos=pos;
	}
	/** Gets the last waiting point's position relative to the sign*/
	public int getPrevWaitPos() { return prevWaitPos; }
	/** Gets the last waiting point's sign*/
	public int getPrevWaitTl() { return prevWaitTl; }
	/** Gets the last waiting point's sign*/
	public boolean getPrevWaitTlColor() { return prevWaitTlColor; }
	/** Gets the current waiting point's sign*/
	public int getCurrentWaitTl() { return waitTl; }
	/** Gets the current waiting point's position relative to the sign*/
	public int getCurrentWaitPos() { return waitPos; }
	/** Gets the current waiting point's sign*/
	public boolean getCurrentWaitTlColor() { return waitTlColor; }

	public boolean didVote() { return voted; }
	public void setVoted(boolean v) { voted = v; }

	public boolean didMove(int cycleNow) {
		cycleAsked = cycleNow;
		return cycleMoved == cycleNow;
	}

	public abstract String getName();
	public abstract int getNumPassengers();
	public abstract int getLength();
	public abstract int getSpeed();
	public abstract int getType();
	public String getVehicleName() { return "unknown"; }
	public String getDriverName() { return "unknown"; }
	public String getDescription() { return "no description available."; }
	public String getPicture() { return null; }
	public String getSound() { return null; }

















	/*============================================*/
	/* Selectable                                 */
	/*============================================*/

	public Rectangle getBounds() { return new Rectangle(0, 0, 0, 0); }
	public Shape getComplexBounds() { return getBounds(); }
	public int getDistance(Point p) { return 0; }
	public Point getSelectionPoint() { return new Point(0, 0); }
	public Point[] getSelectionPoints() { return new Point[0]; }
	public Point getCenterPoint() { return new Point(0, 0); }
	public boolean isSelectable() { return false; }
	public boolean hasChildren() { return false; }
	public Enumeration getChildren() { return null; }
	
	public void paint(Graphics g) throws GLDException
	{
		paint(g, 0, 0, 1.0f);
	}
	
	public void paint(Graphics g, int dx, int dy, double angle) throws GLDException {
		paint(g, dx, dy, 1.0f, angle);
	}
	public abstract void paint(Graphics g, int dx, int dy, float zf) throws GLDException;
	public abstract void paint(Graphics g, int dx, int dy, float zf, double angle) throws GLDException;



















	/*============================================*/
	/* Load/save                                  */
	/*============================================*/
	
 	public void load (XMLElement myElement,XMLLoader loader) throws XMLTreeException,IOException,XMLInvalidInputException
 	{ 	loadData.startNodeId=myElement.getAttribute("start-node-id").getIntValue();
   		loadData.destNodeId=myElement.getAttribute("dest-node-id").getIntValue();
   		cycleMoved=myElement.getAttribute("cycle-moved").getIntValue();
   		cycleAsked=myElement.getAttribute("cycle-asked").getIntValue();
   		position=myElement.getAttribute("position").getIntValue();
   		drivelaneStartTime=myElement.getAttribute("lane-start-time").getIntValue();
   		delay=myElement.getAttribute("delay").getIntValue();
   		prevSign=myElement.getAttribute("prevSign").getIntValue();
			color=new Color (myElement.getAttribute("color-red").getIntValue(),
								myElement.getAttribute("color-green").getIntValue(),
								myElement.getAttribute("color-blue").getIntValue());			
	}	

 	public XMLElement saveSelf () throws XMLCannotSaveException
 	{ 	XMLElement result=new XMLElement(XMLUtils.getLastName(getXMLName()));
   		result.addAttribute(new XMLAttribute("start-node-id",startNode.getId()));
   		result.addAttribute(new XMLAttribute("dest-node-id",destNode.getId()));
   		result.addAttribute(new XMLAttribute("cycle-moved",cycleMoved));
   		result.addAttribute(new XMLAttribute("cycle-asked",cycleAsked));
   		result.addAttribute(new XMLAttribute("position",position));
   		result.addAttribute(new XMLAttribute("type",getType())); // For tunnels
   		result.addAttribute(new XMLAttribute("lane-start-time",drivelaneStartTime));
   		result.addAttribute(new XMLAttribute("delay",delay));
   		result.addAttribute(new XMLAttribute("prevSign",prevSign));
		result.addAttribute(new XMLAttribute("color-green",color.getGreen()));
		result.addAttribute(new XMLAttribute("color-blue",color.getBlue()));
		result.addAttribute(new XMLAttribute("color-red",color.getRed()));
   		return result;
 	}
  
 	public void saveChilds (XMLSaver saver) throws XMLTreeException,IOException,XMLCannotSaveException
 	{ 	// Roadusers don't have child objects
 	}

 	public void setParentName (String parentName)
 	{ 	this.parentName=parentName;
 	}
 
 	class TwoStageLoaderData
 	{ 	int startNodeId,destNodeId;
 	}
 
 	public void loadSecondStage (Dictionary dictionaries) throws XMLInvalidInputException,XMLTreeException
 	{	Dictionary nodeDictionary=(Dictionary)(dictionaries.get("node"));
   		startNode=(Node)(nodeDictionary.get(new Integer(loadData.startNodeId)));
   		destNode=(Node)(nodeDictionary.get(new Integer(loadData.destNodeId)));
 	}

}