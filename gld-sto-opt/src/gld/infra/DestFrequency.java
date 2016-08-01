
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

import gld.infra.*;
import gld.xml.*;
import java.io.IOException;

/**
* Contains a destination frequency for a certain roaduser type.
*
* @author Group Datastructures
* @version 1.0
*/
public class DestFrequency implements XMLSerializable
{
	public int ruType;
	public float freq;
	protected String parentName="model.infrastructure.node";

	/** Empty constructor for loading */	
	public DestFrequency ()
	{ // For loading
	}
	
	/** 
	* Creates an instance initiated with given parameters.
	* @param _ruType Roaduser type.
	* @param _freq Initial frequency.
	*/
	public DestFrequency(int _ruType, float _freq) 
	{	ruType = _ruType;
		freq = _freq;			
	}
	
	// XML Serializable implementation
	
	public void load (XMLElement myElement,XMLLoader loader) throws XMLTreeException,IOException,XMLInvalidInputException
	{ 	ruType=myElement.getAttribute("ru-type").getIntValue();
	  	freq=myElement.getAttribute("freq").getFloatValue();
	}
	public XMLElement saveSelf () throws XMLCannotSaveException
	{ 	XMLElement result=new XMLElement("destfreq");
	  	result.addAttribute(new XMLAttribute("ru-type",ruType));
	  	result.addAttribute(new XMLAttribute("freq",freq));
	 	return result;
	}
 
	public void saveChilds (XMLSaver saver) throws XMLTreeException,IOException,XMLCannotSaveException
	{ 	// A destfrequency has no child objects
	}
 
	public String getXMLName ()
	{ 	return parentName+".destfreq";
	}
	
	public void setParentName (String newParentName)
	{	this.parentName=parentName; 
	}
	
}