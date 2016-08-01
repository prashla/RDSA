
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

import gld.xml.*;
import java.io.IOException;
import java.util.Stack;
import java.util.Vector;

// A very simple number dispenser for ID's. It sucks. It also works.

public class NumberDispenser implements XMLSerializable
{ 	Stack stack;
  	int counter;
	String parentName="model";
  
  	public NumberDispenser ()
  	{ 	stack=new Stack();
    		counter=0;
  	}
  
  	public int get ()
  	{ 	if (stack.isEmpty())
    			return counter++;
    		else
    			return ((Integer)(stack.pop())).intValue();
	}
  
  	public void giveBack (int number)
  	{ 	if (number == counter - 1)
    			counter--;
    		else if (number < counter && ! stack.contains(new Integer (number)))
    			stack.push (new Integer (number));
  	}
	
	public void load (XMLElement myElement,XMLLoader loader) throws XMLTreeException,IOException,XMLInvalidInputException
	{ 	stack=new Stack();
		stack.addAll((Vector)(XMLArray.loadArray(this,loader)));
		counter=myElement.getAttribute("counter").getIntValue();
	}

	public XMLElement saveSelf () throws XMLCannotSaveException
	{ 	XMLElement result=new XMLElement("dispenser");
		result.addAttribute(new XMLAttribute ("counter",counter));
	  	return result;
	}
  
	public void saveChilds (XMLSaver saver) throws XMLTreeException,IOException,XMLCannotSaveException
	{ 	XMLArray.saveArray(stack,this,saver,"stack");
	}

	public String getXMLName ()
	{ 	return parentName+".dispenser";
	}
	
	public void setParentName (String parentName)
	{	this.parentName=parentName; 
	}
	
}
