
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

import java.util.*;

/** This class is an Enumeration which works on objects that implement
  * java.util.List
 */
 
public class ListEnumeration implements Enumeration
{ private Iterator i;

  /** Make a new ListEnumeration
    * @param list The list to use
   */
  public ListEnumeration (List list)
  { i=list.iterator();
  }
  
  // Enumeration implementation
  public boolean hasMoreElements ()
  { return i.hasNext();
  }
  
  public Object nextElement ()
  { return i.next();
  }

}