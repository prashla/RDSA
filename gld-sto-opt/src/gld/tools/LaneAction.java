
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
 * This implements the add lane user action
 *
 * @author Group GUI
 * @version 1.0
 */


public class LaneAction implements ToolAction
{
	protected EditModel model;

	protected Node alphaNode = null;
	protected Point startPoint = null;
	protected Point mousePoint = null;

	public LaneAction(EditModel em) {
		model = em;
	}

	public boolean beingUsed() { return alphaNode != null; }
	
	public void reset() {
		alphaNode = null;
		startPoint = null;
		mousePoint = null;
	}
	
	/**
	 * Starts the action by finding the alpha node of the road
	 *
	 * @param view The main view. Needed to convert coordinates
	 * @param p The point to search for the first node
	 * @return False if no node was found at the given point
	 *         (and the action could not be started consequently)
	 */
	public boolean startAction(View view, Point p) {
		Class[] sf = { Node.class };
		Node clicked = (Node)Selection.selectObject(sf, model.getInfrastructure(), p);
		
		if (clicked != null) {
			alphaNode = clicked;
			startPoint = alphaNode.getCoord();
			mousePoint = p;
			return true;
		}
		return false;
	}
	
	/**
	 * Looks for the second node at given point. If it is found,
	 * ends the action by adding a drivelane to the road connecting the
	 * first and second node.
	 *
	 * @param p The point
	 * @return False if the action was ended without adding a drivelane.
	 */
	public boolean endAction(View view, Point p) {
		Class[] sf = { Node.class };
		Node clicked = (Node)Selection.selectObject(sf, model.getInfrastructure(), p);

		try {
		
			if (clicked == alphaNode) return true; // clicking the alpha node again does not end the action
			if (clicked == null) { reset(); return false; }
			else {
				
				// add drivelane
				Node betaNode = clicked;
				Road[] alphaRoads = alphaNode.getAllRoads();
				Road[] betaRoads = betaNode.getAllRoads();
				Road road = null;
				for (int i=0; i < alphaRoads.length; i++)
					if (alphaRoads[i] != null)
						for (int j=0; j < betaRoads.length; j++)
							if (alphaRoads[i] == betaRoads[j]) road = alphaRoads[i];
				
				if (road == null || road.getNumInboundLanes(betaNode) >= 4) { reset(); return false; }
				
				model.addLane(new Drivelane(road), road, betaNode);

				reset();
			}
		}
		catch (InfraException e) { reset(); Controller.reportError(e); return false; }
		return true;
	}
	
	
	/**
	 * Moves the current mouse point
	 *
	 * @param p The new position of the mouse cursor
	 */
	public void moveAction(View view, Point p) {
		mousePoint = p;
	}
	
	/**
	 * Paints a graphical representation of the status of this action
	 *
	 * @param g The Graphics object to paint on
	 */
	public void paint(Graphics g) {
		g.setXORMode(Color.darkGray);
		g.setColor(Color.lightGray);
		g.drawLine(startPoint.x, startPoint.y, mousePoint.x, mousePoint.y);
	}
}