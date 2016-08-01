
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

package gld.algo.tlc;

import gld.*;
import gld.algo.tlc.*;
import gld.infra.*;
import gld.xml.*;
import java.io.IOException;
import java.util.Dictionary;

/**
 *
 * This class holds a tuple of a TrafficLight and a float value to represent 
 * the reward (Q) for the TrafficLight to be kept.
 *
 * Note : TLControllers typically don't use the XMLSerializable interface
 *        of this class, because they can load it faster themselves.
 * @author Group Algorithms
 * @version 1.0
 */
public class TLDecision implements XMLSerializable,TwoStageLoader
{	private TrafficLight trafficLight;
	private float gain;
	
	private int numGreenCycles=0;
	private int numRedCycles=0;
	private int kappa=0;
	
	protected TwoStageLoaderData loadData=new TwoStageLoaderData();
	protected String parentName="model.tlc";
	
	/** Empty constructor for loading */
	public TLDecision ()
	{}
	
	/**
	 * The constructor for TLDecision.
	 * @param tl The Trafficlight.
	 * The reward value (Q value).
	 */
	public TLDecision(TrafficLight tl, float f)
	{	trafficLight = tl;
		gain = f;
	}
	
	/**
	 * Returns the TrafficLight.
	 * @return The TrafficLight.
	 * @see gld.infra.TrafficLight
	 */	
	public TrafficLight getTL() {
		return trafficLight;
	}
	
	/**
	 * Returns the Q value.
	 * @return The qValue.
	 */	
	public float getGain() {
		return gain;
	}
	
	public void setGain(float _g) {
		gain = _g;
	}
	
	public int getNumGreenCycles() { return numGreenCycles; }
	public void setNumGreenCycles(int k) { numGreenCycles = k;}
	public void addNumGreenCycles(int k) {numGreenCycles += k;}
	public int getNumRedCycles() { return numRedCycles; }
	public void setNumRedCycles(int k) { numRedCycles = k;}
	public void addNumRedCycles(int k) {numRedCycles += k;}
	
	//For sotl request
	public int getKappa() { return kappa; }
	public void setKappa(int k) { kappa = k;}
	public void addKappa(int k) {kappa += k;}
	
	public void load (XMLElement myElement,XMLLoader loader) throws XMLTreeException,IOException,XMLInvalidInputException 
	{	gain=myElement.getAttribute("q").getFloatValue();
		loadData.tlId=myElement.getAttribute("tl-id").getIntValue();
	}
		      
	public XMLElement saveSelf () throws XMLCannotSaveException { 	
		XMLElement result=new XMLElement("decision");		
		result.addAttribute(new XMLAttribute("q",gain));
		result.addAttribute(new XMLAttribute("tl-id",trafficLight.getId()));
	  	return result;
	}
  
	public void saveChilds (XMLSaver saver) throws XMLTreeException,IOException,XMLCannotSaveException 
	{	// A TLDecision has no child objects
	}

	public String getXMLName () {
		return parentName+".decision";
	}
	
	public void setParentName (String parentName)
	{	this.parentName=parentName;
	}
	
	class TwoStageLoaderData 
	{	int tlId;
	}
		
	public void loadSecondStage (Dictionary dictionaries) throws XMLInvalidInputException,XMLTreeException {
		trafficLight=(TrafficLight)((Dictionary)(dictionaries.get("sign"))).get
						(new Integer(loadData.tlId));		
	}
}
