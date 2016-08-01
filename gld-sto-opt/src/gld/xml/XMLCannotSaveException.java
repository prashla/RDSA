
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

import gld.GLDException;

/** Thrown by a XMLSerializable when it is told to save it self. But that
  * cannot be done. (e.g. because the object is in an inconsistent state)
 */
 
public class XMLCannotSaveException extends GLDException
{ /** Oops... 
    * @param message A description of the problem.
   */
  public XMLCannotSaveException (String message)
  { super(message);
  }
}