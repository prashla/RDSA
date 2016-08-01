
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
 * This implements scroll user action (like in acrobat reader)
 *
 * @author Group GUI
 * @version 1.0
 */


public class ScrollAction implements ToolAction
{
	protected Point startPoint = null;
	protected Point startScrollPoint = null;
	protected ViewScroller viewScroller = null;
	protected boolean skipnext = false;

	public ScrollAction(ViewScroller vs) {
		viewScroller = vs;
	}

	public boolean beingUsed() { return startPoint != null; }

	public boolean startScroll(View view, Point p)
	{
		startPoint = p;
		startScrollPoint = viewScroller.getScrollPosition();
		return true;
	}

	public void doScroll(View view, Point p)
	{
		if (startPoint.distance(p) > 5)
		{
			if (skipnext) { skipnext = false; return; }
			int offx = p.x - startPoint.x;
			int offy = p.y - startPoint.y;
			Point to = new Point(startScrollPoint.x - offx, startScrollPoint.y - offy);
			startScrollPoint = to;
			viewScroller.setScrollPosition(to);
			skipnext = true;
		}
	}

	public void endScroll(View view, Point p)
	{
		Point to = new Point(startScrollPoint.x + p.x - startPoint.x, startScrollPoint.y + p.y - startPoint.y);
		viewScroller.setScrollPosition(to);
		startPoint = null;
		startScrollPoint = null;
	}
}