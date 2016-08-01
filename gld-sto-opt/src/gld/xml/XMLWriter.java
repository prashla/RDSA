
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

import gld.utils.StringUtils;
import java.io.*;

/** A utility class to write XML to an OutputStream
 */


class XMLWriter extends PrintWriter
{ 	OutputStream stream;

  /** Make a new XMLWriter
   * @param stream The OutputStream on which this XMLWriter is based
	* @param autoFlush Indicates if the internal PrintWriter needs autoflush
   *
  */
  public XMLWriter (OutputStream stream,boolean autoFlush)
  { 	super(stream,autoFlush);
  		this.stream=stream;
  }
  
  /** Make a new XMLWriter
   * @param stream The OutputStream on which this XMLWriter is based
   *
  */
  public XMLWriter (OutputStream stream)
  { super(stream);
  }
  
  
  /** Write the open tag of a XMLElement 
    * @param element The XMLElement
	 * @param indent The number of spaces to indent the tag
	*/
  public void writeOpenTag (XMLElement element,int indent)
  { 	println(StringUtils.repeat(' ',indent)+element.getOpenTag());
  }
	
  /** Write the open tag of a XMLElement 
    * @param element The XMLElement
	*/
  public void writeOpenTag (XMLElement element)
  {	writeOpenTag(element,0);
  }
  
  /** Write the close tag of a XMLElement 
    * @param element The XMLElement
	 * @param indent The number of spaces to indent the tag
	*/
  public void writeCloseTag (XMLElement element,int indent)
  { 	println(StringUtils.repeat(' ',indent)+element.getCloseTag());
  }
	
  /** Write the close tag of a XMLElement 
    * @param element The XMLElement
	*/
  public void writeCloseTag (XMLElement element)
  {	writeCloseTag(element,0);
  }
  
  /** Write both open and close tags of a XMLElement 
    * @param element The XMLElement
	 * @param indent The number of spaces to indent the tags
	*/
  public void writeAtomaryElement (XMLElement element,int indent)
  { 	println(StringUtils.repeat(' ',indent)+
      	element.getOpenTag()+" "+element.getCloseTag());
  }
	
  /** Write both open and close tags of a XMLElement 
    * @param element The XMLElement
	*/
  public void writeAtomaryElement (XMLElement element)
  {	writeAtomaryElement(element,0);
  }
  
}
