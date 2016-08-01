
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
import java.io.IOException;
import java.util.Dictionary;

/**
 *
 * Sign class.
 *
 * @author Group Datastructures
 * @version 1.0
 */

public abstract class Sign implements XMLSerializable, TwoStageLoader
{
	protected boolean state = false;
	protected Node node;
	protected Drivelane lane;
	/** Data for loading the second stage */
	protected TwoStageLoaderData loadData = new TwoStageLoaderData();
	protected String parentName="model.infrastructure.lane";

	public static final int SIGN = 0;
	public static final int TRAFFICLIGHT = 1;
	public static final int NO_SIGN = 2;
	
	/** last cycle when switched green/red */
	public int cycleSwitchedToGreen =-10;
	public int cycleSwitchedToRed =-10;


	public Sign(Node _node, Drivelane _lane) {
		node = _node;
		lane = _lane;
		state = false;
	}
	
	public Sign () { }
	
	/** Smallest factory ever */
	public static Sign getInstance(int type) throws InfraException {
		if (type == TRAFFICLIGHT) return new TrafficLight();
		if (type == NO_SIGN) return new NoSign();
		throw new InfraException("Unknown sign type: " + type);
	}

	/** This will reset the sign to its default state (false) */
	public void reset() {
		state = false;
	}

	/** Returns true if this Sign should be handled by an external algorithm (TC-3 for example) */
	public abstract boolean needsExternalAlgorithm();
	
	/**
	 * Returns true if the Roaduser at the start of the Drivelane may cross the Node.
	 * Default behavior is implemented to return state.
	 * Any Sign using an external algorithm does not need to override this method.
	 */
	public boolean mayDrive() { return state; }

	/** Returns the Id of this sign. It is the same as the Id of the drivelane this sign is on. */
	public int getId() { return lane.getId(); }

	/** Returns the current state */
	public boolean getState() { return state; }
	/** Sets the current state */
	public void setState(boolean b) { state = b; }

	/** Returns the Drivelane this Sign is on */
	public Drivelane getLane() { return lane; }
	/** Sets the Drivelane this Sign is on */
	public void setLane(Drivelane l) { lane = l; }
	
	/** Returns the Node this Sign is on */
	public Node getNode() { return node; }
	/** Sets the Node this Sign is on */
	public void setNode(Node n) { node = n; }

	/** Returns the type of this Sign */
	public abstract int getType();
	
	/** returns the last cycle when sign was switched green */
	public int getCycleSwitchedToGreen() { return cycleSwitchedToGreen; }
	/** sets the cycle when sign is switched green */
	public void setCycleSwitchedToGreen(int c) { cycleSwitchedToGreen = c; }
	
	/** returns the last cycle when sign was switched red */
	public int getCycleSwitchedToRed() { return cycleSwitchedToRed; }
	/** sets the cycle when sign is switched green */
	public void setCycleSwitchedToRed(int c) { cycleSwitchedToRed = c; }
	
	// Common XMLSerializable implementation for subclasses
	
	public void load (XMLElement myElement,XMLLoader loader) throws XMLTreeException,IOException,XMLInvalidInputException
	{
		state=myElement.getAttribute("state").getBoolValue();
		loadData.nodeId=myElement.getAttribute("node-id").getIntValue();
	}

	public XMLElement saveSelf () throws XMLCannotSaveException
	{ 	XMLElement result=new XMLElement("sign");
	  	result.addAttribute(new XMLAttribute("type",getType()));
	  	result.addAttribute(new XMLAttribute("state",state));
	  	result.addAttribute(new XMLAttribute("node-id",node.getId()));
	  	// Lane id doesn't have to be saved, because it is set by the parent
	  	// lane on loading
	 	return result;
	}
  
	public void saveChilds (XMLSaver saver) throws XMLTreeException,IOException,XMLCannotSaveException
	{
	}

	public void setParentName (String parentName)
	{	this.parentName=parentName; 
	}

	class TwoStageLoaderData 
 	{ 	int nodeId;
 	}
 
 	public void loadSecondStage (Dictionary dictionaries) throws XMLInvalidInputException,XMLTreeException
 	{ 	node=(Node)(((Dictionary)(dictionaries.get("node"))).get(new Integer(loadData.nodeId)));
 		//System.out.println("Node gotten:"+node);
 	}
 }
