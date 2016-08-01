
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

package gld.edit;

import gld.*;

import java.awt.*;

/**
*
* Overlay for <code>gld.View</code>. Shows a grid.
*
* @author Group GUI
* @version 1.0
*/

public class GridOverlay implements Overlay
{
	protected final static int GRID_SPACING = 20;

	Dimension gridSize;	
	
	/** Creates a default <code>GridOverlay</code>. */
	public GridOverlay(Dimension size)
	{
		gridSize = size;
	}
	
	public int overlayType() { return 2; }
	
	/** Paints the grid. */
	public void paint(Graphics g)
	{
		g.setPaintMode();
		g.setColor(new Color(128, 128, 128, 128));
		
		int top = -(int)(gridSize.height / 2);
		int left = -(int)(gridSize.width / 2);
		top -= top % GRID_SPACING;
		left -= left % GRID_SPACING;
		
		int cellSize = GRID_SPACING;
		
		for(int x = 0; x < gridSize.width; x += cellSize)
			g.drawLine(x + left, top, x + left, top + gridSize.height);
		
		for(int y = 0; y < gridSize.height; y += cellSize)
			g.drawLine(left, y + top, left + gridSize.width, y + top);
	}
}