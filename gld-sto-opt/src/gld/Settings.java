
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

package gld;

import gld.xml.*;
import java.io.IOException;
import java.util.NoSuchElementException;

/** General interface for classes which hold settings */

public interface Settings extends XMLSerializable
{	
	/** Load settings from settings file */
	public void doLoad () throws IOException,XMLTreeException,XMLInvalidInputException;
	
	/** Save settings to settings file */
	public void doSave () throws IOException,XMLTreeException,XMLCannotSaveException;
	
	/** Indicates if this Settings are predefined. If so, then the user cannot
	  * add new properties to this object. So, in that case, it is not possible
	  * to do a setProperty with a name that does not exist yet.
	  * So it is OK to change existing properties of a predefined settings object,
	  * but you cannot add new properties.
	 */
	public boolean isPredefined ();
	
	/** Set a string property to a certain value 
	  * @param name The name of the property
	  * @param value The new value of the property
	  * @param NoSuchElementException If the settings are predefined and there
	  *        is no property with the specified name
	 */
	public void setProperty(String name,String value) throws	NoSuchElementException;

	/** Set an int property to a certain value 
	  * @param name The name of the property
	  * @param value The new value of the property
	  * @param NoSuchElementException If the settings are predefined and there
	  *        is no property with the specified name 
	 */
	public void setProperty(String name,int value) throws	NoSuchElementException;	
	
	/** Set a boolean property to a certain value 
	  * @param name The name of the property
	  * @param value The new value of the property
	  * @param NoSuchElementException If the settings are predefined and there
	  *        is no property with the specified name
	 */
	public void setProperty(String name,boolean value) throws NoSuchElementException;	
	
	/** Gets the String value of a certain property
	  * @param name The name of the property 
	  * @throws NoSuchElementException if there is no property with that name
	 */	
	public String getPropertyStringValue (String name) throws NoSuchElementException;
	
	/** Gets the value of an int property 
	  * @param name The name of the property 
	  * @throws NoSuchElementException if there is no property with that name
	 */	
	public int getPropertyIntValue (String name) throws NoSuchElementException;
	
	/** Gets the value of a boolean property
	  * @param name The name of the property 
	  * @throws NoSuchElementException if there is no property with that name
	 */	
	public boolean getPropertyBooleanValue (String name) throws	NoSuchElementException;
	
	/** Gets the value of a float property
	  * @param name The name of the property 
	  * @throws NoSuchElementException if there is no property with that name
	 */	
	public float getPropertyFloatValue (String name) throws NoSuchElementException;
	
}
