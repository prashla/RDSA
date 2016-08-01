
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

import java.util.Dictionary;

/** This interface can be implemented by XMLSerializable objects if they
  * need additional info from their parent for initializing themselves after
  * they have been loaded by the parser.
 */

public interface TwoStageLoader
  /** This method is called by the parent of the object to provide it with
    * additional information after the XML parser has loaded it via the
    * XMLSerializable interface.
    * @param dictionaries The info is passed in the form of a dictionary
    *                     of dictionaries. Each dictionary is a translation
    * 			  table from an unique id to an object.		 		  
   */ 
    
{ public void loadSecondStage (Dictionary dictionaries) throws XMLInvalidInputException,XMLTreeException;
}
