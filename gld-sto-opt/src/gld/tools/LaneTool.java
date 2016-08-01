
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
import gld.edit.*;

import java.awt.*;


/**
 * This tool implements the LaneAction
 *
 * @author Group GUI
 * @version 1.0
 */

public class LaneTool extends PopupMenuTool
{
	protected LaneAction la;
	
	public LaneTool(EditController c) {
		super(c);
		la = new LaneAction(c.getEditModel());
	}

	public void mousePressed(View view, Point p, Tool.Mask mask)
	{
		if (la.beingUsed()) {
			if (mask.isRight()) {
				la.reset();
				view.repaint();
			}
			return;
		}

		super.mousePressed(view, p, mask);

		if (!mask.isLeft()) return;
		if (la.startAction(view, p)) view.repaint();
	}

	public void mouseReleased(View view, Point p, Tool.Mask mask)
	{
		if (!mask.isLeft() && !la.beingUsed()) return;
		la.endAction(view, p);
		view.repaint();
	}

	public void mouseMoved(View view, Point p, Tool.Mask mask) {
		if (la.beingUsed()) {
			la.moveAction(view, p);
			view.repaint();
		}
	}

	public int overlayType() { return 1; }
	
	public void paint(Graphics g) throws GLDException {
		if (la.beingUsed()) la.paint(g);
	}
	
	public Panel getPanel() { return new Panel(null); }
}