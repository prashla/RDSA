
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
 * You click on two Nodes with a RoadTool to create a Road.
 *
 * @author Group GUI
 * @version 1.0
 */

public class RoadTool extends PopupMenuTool
{
	protected RoadAction ra;
	
	public RoadTool(EditController c) {
		super(c);
		ra = new RoadAction(c.getEditModel());
	}

	public void mousePressed(View view, Point p, Tool.Mask mask)
	{
		if (ra.beingUsed()) {
			if (mask.isRight()) {
				ra.reset();
				view.repaint();
			}
			return;
		}

		super.mousePressed(view, p, mask);

		if (!mask.isLeft()) return;
		if (ra.startAction(view, p)) view.repaint();
	}

	public void mouseReleased(View view, Point p, Tool.Mask mask)
	{
		if (!mask.isLeft() && !ra.beingUsed()) return;
		ra.nextAction(view, p);
		view.repaint();
	}

	public void mouseMoved(View view, Point p, Tool.Mask mask) {
		if (ra.beingUsed()) {
			ra.moveAction(view, p);
			view.repaint();
		}
	}

	public int overlayType() { return 1; }
	
	public void paint(Graphics g) throws GLDException {
		if (ra.beingUsed()) ra.paint(g);
	}
	
	public Panel getPanel() { return new Panel(null); }
}