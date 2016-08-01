
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
import java.util.*;

/**
 *
 * Sign implementing basic trafic rules.
 *
 * @author Group Datastructures
 * @version 1.0
 */

public class NoSign extends Sign
{
	protected final static int type = Sign.NO_SIGN;
    
	public NoSign(Node _node, Drivelane _lane ) {
		super(_node, _lane);
	}
	
	public NoSign () { }	

	/** Returns the type of this Sign */
	public int getType() { return type; }

	public boolean needsExternalAlgorithm() { return false; }
	
	public boolean mayDrive() { return true; }

	
	// Specific XMLSerializable implementation 

	public String getXMLName ()
	{ 	return parentName+".sign-no";
	}
	
	
	public XMLElement saveSelf () throws XMLCannotSaveException
	{ 	XMLElement result=super.saveSelf();
		result.setName("sign-no");
		return result;
	}
	
}
