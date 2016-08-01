
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

import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.util.*;

import gld.infra.*;
import gld.tools.Tool;
import gld.xml.*;

/**
 *
 * This view paints an entire infrastructure
 *
 * @author  Group Interfaces
 * @version 1.0
 */

public class MainView extends View
	
{
	/** The current infrastructure */
	protected Infrastructure infra;
	
	/**
	* Creates a new MainView
	* @param i The Infrastructure to be displayed
		*/
	public MainView(Infrastructure i) {
		super(i.getSize());
		infra = i;
	}

	/** Returns the infrastructure 
	* @return infra The infrastructure
	*/
	public Infrastructure getInfrastructure() { return infra; }
	/** Sets the infrastructure */
	public void setInfrastructure(Infrastructure i) { infra = i; }


	/** Draw a buffer for the current infrastructure*/
	public void fillBuffer(Graphics2D g) {
		Node[] nodes = infra.getAllNodes();
		Road[] roads;

		try {
			for (int i=0; i < nodes.length; i++) {
				nodes[i].paint(g);
				roads = nodes[i].getAlphaRoads();
				for (int j=0; j < roads.length; j++) {
					try {
						roads[j].paint(g);
					}
					catch (ConcurrentModificationException e) {
						try {
							roads[j].paint(g); // try once more
						}
						catch (ConcurrentModificationException ex) { }
					}
				}
			}
		}
		catch (GLDException e) {
			Controller.reportError(e);
		}
	}
}