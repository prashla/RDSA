
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
 * Left-click to zoom in, right-click to zoom out. That is, when this <code>Tool</code>
 * is the currently selected <code>Tool</code>.
 *
 * @author Group GUI
 * @version 1.0
 */

public class ZoomTool implements Tool
{
	protected ZoomAction za;

	/**
	 * Creates a <code>ZoomTool</code>.
	 *
	 * @param c The <code>Controller</code> controlling this <code>Tool</code>.
	 */
	public ZoomTool(Controller c) {
		za = new ZoomAction(c);
	}

	/**
	 * Invoked when the user releases a mouse button.
	 * Zoom in on left-click, zoom out on right-click.
	 *
	 * @param view The <code>View</code> that the event originates from.
	 * @param p The coordinates in the infrastructure the mouse cursor was at when the event was generated.
	 * @param mask Identifies which button was pressed, as well as any aditional sytem keys
	 */
	public void mousePressed(View view, Point p, Tool.Mask mask)
	{
		if (mask.isLeft()) za.doZoom(view, p, ZoomAction.IN);
		else if (mask.isRight()) za.doZoom(view, p, ZoomAction.OUT);
	}

	public void mouseReleased(View view, Point p, Tool.Mask mask) { }
	public void mouseMoved(View view, Point p, Tool.Mask mask) { }
	public int overlayType() { return 0; }
	public void paint(Graphics g) throws GLDException { }

	public Panel getPanel() { return new Panel(null); }
}