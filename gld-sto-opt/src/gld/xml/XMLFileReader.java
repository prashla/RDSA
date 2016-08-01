
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

import java.util.zip.GZIPInputStream;
import java.io.*;

/** A XMLReader which can read files
 */

class XMLFileReader extends XMLReader
{  /** Make a new XMLFileReader
    * @param file The file to read the data from
    * @throws IOException If the file cannot be opened because of an IO error
   */
  public XMLFileReader (File file) throws IOException
  { // Just to gain some speeeeed! dropped the 'new GZIPInputStream('
	super(new FileInputStream(file));
  }

}