
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


/** The GLD user documentation defines InstantiationAssistant as 
  * "Your friendly class creator who gives you a warm and fuzzy feeling
  * inside". The technical documentation says that it is just another 
  * interface which classes can implement if they can create classes that
  * the XML Parser cannot instantiate. For instance non-static inner classes.
  * The parser has to figure out which InstantiationAssistant to use : for an
  * inner class it will normally use the parent class.
 */

public interface InstantiationAssistant
{	/** Tell our parser if we can create an instance of a certain class
	  * @param request The class of the object that the parser wants
	  * @returns Can we create that?
	 */
	public boolean canCreateInstance (Class request);

	/** Create an instance of a certain class for our parser
	  * @param request The class of the object that the parser wants
	  * @returns An instance of that object
	  * @throws ClassNotFoundException If we don't know that class
	  * @throws InstantiationException In case something goes wrong
	 */
	public Object createInstance (Class request) 
	throws ClassNotFoundException,InstantiationException,IllegalAccessException;
	
}
