
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
 * Tool for selecting through clicks or drags Nodes and Roads.
 *
 * @author Group GUI
 * @version 1.0
 */
 
public class SelectTool extends PopupMenuTool
{
	protected SelectAction selectAction;
	
	public SelectTool(Controller con) {
		super(con);
		selectAction = new SelectAction(con.getCurrentSelection());
	}
	
	public void mousePressed(View view, Point p, Tool.Mask mask) {
		super.mousePressed(view, p, mask);
		if (mask.isLeft()) selectAction.doStart(view, p);
	}

	public void mouseReleased(View view, Point p, Tool.Mask mask) {
		if (selectAction.beingUsed()) {
			int type = SelectAction.NEW;
			if (mask.isControlDown()) type = SelectAction.INVERT;
			if (mask.isShiftDown()) type = SelectAction.ADD;
			selectAction.endDrag(view, p, type);
		}
	}

	public void mouseMoved(View view, Point p, Tool.Mask mask) {
		if (selectAction.beingUsed()) selectAction.doDrag(view, p);
	}
	
	public int overlayType() { return 1; }
	
	public void paint(Graphics g) throws GLDException {
		if (selectAction.beingUsed()) selectAction.paint(g);
	}
	
	public Panel getPanel() { return new Panel(null); }
}