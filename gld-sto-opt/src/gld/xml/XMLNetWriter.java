
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

import java.io.IOException;
import java.io.OutputStream;
import java.net.Socket;

/** A XMLWriter which can write to sockets
 */

class XMLNetWriter extends XMLWriter
{  /** Make a new XMLNetWriter
    * @param host The hostname to connect to
	 * @param port The port number to writer the data to
    * @throws IOException Some I/O error
	 * @throws UnknownHostException If the hostname cannot be resolved
   */
  protected XMLNetWriter(String host,int port) throws IOException
  { super((new Socket(host,port)).getOutputStream(),true);
  }
  
}
