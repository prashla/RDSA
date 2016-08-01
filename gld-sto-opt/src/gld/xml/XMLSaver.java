
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

package gld.xml;

import gld.*;
import gld.xml.*;
import gld.utils.*;
import java.io.*;
import java.util.*;

// TODO : - write Junit test

/** This class takes care of saving the state of the program to a XML
  * file. It also provides methods which XMLSerializable objects can
  * use to save their child objects.
 */

public class XMLSaver
{ protected XMLWriter writer;
  protected XMLStack stack;
  
  /** Make a new XMLSaver
   */
  protected XMLSaver ()
  { stack=new XMLStack();
  }

  /** Makes a new XMLSaver which saves to a file
    * @param file The file to save to
    * @throws IOException If the XMLSaver cannot open the XML file because
    *         of an IO error.
   */
  	public XMLSaver (File file) throws IOException
  	{ 	this();
  		writer=new XMLFileWriter(file);
  	}

  /** Makes a new XMLSaver which saves to a socket
    * @param hostname The hostname to connect to
	 * @param port The port to connect to
    * @throws IOException If the XMLSaver cannot open the socket connection
    *         because of an IO error.
   */
 	public XMLSaver (String hostname,int port) throws IOException
  	{ 	this();
		writer=new XMLNetWriter(hostname,port);
  	}
  
  /** A crude indication if we can write to this XMLSaver */
  public boolean hasStream ()
  { return writer != null ;
  }
  
  /** Close this XMLSaver */
  public void close () 
  { 	writer.flush();
  		writer.close();
  		writer=null;
  }

  /** Save a XMLSerializable object to the XML file
    * @param object The object to save
    * @throws IOException If the XMLSaver cannot write to the XML file because
    *         of an IO error.
    * @throws XMLTreeException If something goes wrong while the parser is 
    *         building the XML tree structure for this object
    * @throws XMLCannotSaveException If this object cannot save itself
   */
  public void saveObject (XMLSerializable object) throws IOException,XMLTreeException,XMLCannotSaveException
  { //Initialize various names of the object
    String fullName  =object.getXMLName(),
           parentName=XMLUtils.getParentName(object),
           firstName =XMLUtils.getLastName(object);
    // Check if the object can be saved.
    stack.assertIsCurrentParent(parentName);
    // Save element of object
    stack.push(firstName);
    XMLElement element=object.saveSelf();
	 writer.writeOpenTag(element,indent());
    // Let object save its child objects
    object.saveChilds(this);
    // Finish element
	 writer.writeCloseTag(element,indent());
    stack.pop();
    // Check if the stack is OK
    if ( (! stack.getBranchName().equals(parentName)) &&
	 (! XMLUtils.getGenericName(stack.getBranchName()).equals
           (XMLUtils.getGenericName(parentName))))	    
    { throw new XMLTreeException
      ("The XMLSerializable "+fullName+
       " (or possibly one of its child objects) left the XML tree "+
       "in an invalid state. Terminating save.\n"+
       "Expected branchName :"+stack.getBranchName()+"/"+
                               XMLUtils.getGenericName(stack.getBranchName())+"\n"+
       "but branchname was  :"+parentName+"/"+
                               XMLUtils.getGenericName(parentName));
    }
  }

 /** Save an enum_vareration of XML Serializables to file 
    * @param e The enum_vareration of XML Serializables
    * @throws IOException If the XMLSaver cannot write to the XML file because
    *         of an IO error.
    * @throws XMLTreeException If something goes wrong while the parser is 
    *         building the XML tree structure for one of the objects.
    * @throws XMLCannotSaveException If one of the objects cannot save
    *         itself for one reason or another.
    * @throws ClassCastException If the enum_vareration contains an object that
    *         is not XMLSerializable
   */
  public void saveEnumerationObjects (Enumeration e) throws ClassCastException,XMLTreeException,IOException,XMLCannotSaveException
  { while (e.hasMoreElements())
          saveObject((XMLSerializable)(e.nextElement()));
  }
  
 /** Save an atomary XML element
    * @param el The XMLElement to save
    * @throws IOException If the XMLSaver cannot write to the XML file because
    *         of an IO error.
    * @throws XMLTreeException If there is a parser problem with writing
    *         the element.
   */
  public void saveAtomaryElement (XMLSerializable parent,XMLElement el) throws IOException,XMLTreeException
  { 	stack.assertIsCurrentParent(parent);
  		writer.writeAtomaryElement(el,indent());
  }
  
  /** @return An internal indentation string */
  protected int indent ()
  { return stack.size();
  }
  

}
