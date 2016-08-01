
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
 * Uses the ScrollAction to scroll the view
 *
 * @author Group GUI
 * @version 1.0
 */

public class ScrollTool extends PopupMenuTool
{
	protected ScrollAction sa;

	public ScrollTool(Controller c) {
		super(c);
		sa = new ScrollAction(c.getViewScroller());
	}

	public void mousePressed(View view, Point p, Tool.Mask mask)
	{
		super.mousePressed(view, p, mask);
		if (mask.isLeft()) sa.startScroll(view, p);
	}

	public void mouseMoved(View view, Point p, Tool.Mask mask) {
		if (mask.isLeft()) sa.doScroll(view, p);
	}
	public void mouseReleased(View view, Point p, Tool.Mask mask) {
		if (mask.isLeft()) sa.endScroll(view, p);
	}
	public int overlayType() { return 0; }
	public void paint(Graphics g) throws GLDException { }

	public Panel getPanel() { return new Panel(null); }
}