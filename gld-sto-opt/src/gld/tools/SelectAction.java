
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
 * This implements the PopupMenu user action
 *
 * @author Group GUI
 * @version 1.0
 */


public class SelectAction implements ToolAction
{
	public static final int INVERT = 1;
	public static final int ADD = 2;
	public static final int NEW = 3;
	
	Selection currentSelection;
	Point startPoint = null;
	Point endPoint = null;
	Point prevPoint = null;
	boolean isDragging = false;
	
	public SelectAction(Selection s) {
		currentSelection = s;
	}

	/** Returns the enclosing rectangle of the last (current) drag operation */
	private Rectangle getEnclosingRectangle() {
		Rectangle r = new Rectangle(startPoint);
		if (endPoint == null)
			r.add(prevPoint);
		else
			r.add(endPoint);
		return r;
	}

	public boolean beingUsed() { return isDragging; }
	
	public void doStart(View view, Point p) {
		startPoint = p;
		prevPoint = p;
		isDragging = true;
	}
	
	public void doDrag(View view, Point p) {
		prevPoint = endPoint;
		endPoint = p;
		if (startPoint != null) {
			Rectangle r = new Rectangle(view.toView(startPoint));
			if (endPoint != null) r.add(view.toView(endPoint));
			if (prevPoint != null) r.add(view.toView(prevPoint));
			r.grow(1, 1);
			view.repaint(r.x, r.y, r.width, r.height);
		}
	}
	
	public void endDrag(View view, Point p, int type) {
		endPoint = p;
		isDragging = false;
		Rectangle r = getEnclosingRectangle();
		if (r.width < 5 && r.height < 5) {
			Point p2 = new Point(r.x + 2, r.y + 2);
			if (type == INVERT)
				currentSelection.invertWithSelection(p2);
			else if (type == ADD)
				currentSelection.addToSelection(p2);
			else
				currentSelection.newSelection(p2);
		}
		else {
			if (type == ADD || type == INVERT)
				currentSelection.addToSelection(r);
			else
				currentSelection.newSelection(r);
		}
	}

	public void paint(Graphics g) {
		g.setXORMode(Color.darkGray);
		g.setColor(Color.lightGray);
		Rectangle r = getEnclosingRectangle();
		g.drawRect(r.x, r.y, r.width, r.height);
	}
}