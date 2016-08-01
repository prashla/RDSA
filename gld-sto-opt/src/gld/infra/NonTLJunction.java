
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
import gld.utils.*;
import gld.GLDException;
import gld.xml.*;

import java.awt.Point;
import java.awt.Graphics;
import java.awt.Color;
import java.io.IOException;
import java.util.*;

/**
 *
 * Junction without trafficlights, so normal traffic rules are applied
 * @author Group Datastructures
 * @version 1.0
 */


public class NonTLJunction extends Junction
{
	/** The type of this node */
	protected static final int type = Node.NON_TL;
	
	/** A Vector containing all Signs on this node */
	protected NoSign[] signs = { };
	
	/** Creates an empty junction (for loading) */
	 public NonTLJunction () { }
	
	// Guess you want to change this to a new sign?
	public int getDesiredSignType() { return Sign.NO_SIGN; }

	/**
	 * Creates a new standard Junction
	 *
	 * @param _coord The coordinates of this node on the map in pixels.
	 */
	public NonTLJunction(Point _coord) {
		super(_coord);
	}

	// Specific XMLSerializable implementation 

	public String getXMLName ()
	{ 	return parentName+".node-junction-nontl";
	}
	
	public XMLElement saveSelf () throws XMLCannotSaveException
	{ 	XMLElement result=super.saveSelf();
		result.setName("node-junction-nontl");
	  	return result;
	}


}