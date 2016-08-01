
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

import java.io.EOFException;
import java.io.InputStream;
import java.io.IOException;
import java.net.Socket;

/** XMLReader which can read data from a Socket;
 */

class XMLNetReader extends XMLReader
{ Socket socket;
	 /** Make a new XMLNetReader
     * @param socket The socket to read the data from
    */
  public XMLNetReader (Socket socket) throws IOException
  {	super(socket.getInputStream());
  		this.socket=socket;
  }
  
  public void close () throws IOException
  	{	socket.close();
		super.close();
	}
  
}
