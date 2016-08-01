
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

import java.util.zip.GZIPOutputStream;
import java.io.*;

/** A XMLWriter which writes to files
 */

class XMLFileWriter extends XMLWriter
{  /** Make a new XMLFileWriter
    * @param file The file to save the data to
    * @throws IOException If the file cannot be opened because of an IO error
   */
  protected XMLFileWriter(File file) throws IOException
  { 	super(new FileOutputStream(file));
  }

}
