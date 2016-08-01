
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
import gld.infra.*;

import java.awt.*;

/**
 * This implements zoom user action
 *
 * @author Group GUI
 * @version 1.0
 */


public class MoveAction implements ToolAction
{
	protected EditModel model;
	protected Node node = null;         // the Node we are moving
	protected Point startPoint = null;  // in infra space
	protected int xOffset = 0;          // node.coord.x - curpoint.x
	protected int yOffset = 0;          // node.coord.y - curpoint.y

	public MoveAction(EditModel em) {
		model = em;
	}

	public boolean beingUsed() { return node != null; }
	
	public boolean startMove(View view, Point p) {
		Class[] sf = { Node.class };
		Selection s = new Selection(view, sf, model.getInfrastructure());
		s.newSelection(p);
		
		if (!s.isEmpty())
		{
			node = (Node)s.getSelectedObjects().getFirst();
			view.setCursor(Cursor.getPredefinedCursor(Cursor.MOVE_CURSOR));
			startPoint = p;
			xOffset = node.getCoord().x - p.x;
			yOffset = node.getCoord().y - p.y;
			return true;
		}
		return false;
	}
	
	public void doMove(View view, Point p) {
		p.x += xOffset;
		p.y += yOffset;

		if (node != null && p.distance(startPoint) > 4 && nodeMayMove(node, view, p)) {
			model.moveNode(node, p);
		}
	}

	public void endMove(View view, Point p) {
		if (node == null) return;
		
		p.x += xOffset;
		p.y += yOffset;
		
		view.setCursor(null);

		// Check if node is allowed to move here
		if (!nodeMayMove(node, view, p)) {
			node = null;
			return;
		}

		// OK, move node
		model.moveNode(node, p);
		node = null;
		return;
	}

	/** Returns true if given node may be moved to given point */
	public boolean nodeMayMove(Node n, View view, Point p) {
		Class[] sf = { Node.class };
		Selection s = new Selection(view, sf, model.getInfrastructure());
		Rectangle r = new Rectangle(p);
		int growsize = Infrastructure.blockLength * 2 + Infrastructure.blockWidth * 8 + 2;
		r.grow(growsize, growsize);
		s.newSelection(r);
		
		s.deselect(n);

		return s.isEmpty();
	}
}