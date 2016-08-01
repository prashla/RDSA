
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

package gld.tools;

import gld.*;

import java.awt.*;

/**
 * This implements zoom user action
 *
 * @author Group GUI
 * @version 1.0
 */


public class ZoomAction implements ToolAction
{
	public static final int IN = 1;
	public static final int OUT = 2;
	
	Controller controller;
	
	public ZoomAction(Controller con) {
		controller = con;
	}
	
	public boolean beingUsed() { return false; }
	
	public void doZoom(View view, Point p, int type)
	{
		if (type == IN) controller.zoomIn(p);
		else controller.zoomOut(p);
	}
}