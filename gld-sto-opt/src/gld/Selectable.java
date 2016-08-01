
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

package gld;

import gld.GLDException;

import java.awt.*;
import java.util.*;

/**
 *
 *
 * @author Group Datastructures
 * @version 1.0
 */

public interface Selectable extends SelectionStarter
{
	/** Returns the bounding box of this Selectable */
	public Rectangle getBounds();
	/** Returns a complex bounding box of this Selectable */
	public Shape getComplexBounds();
	/** Returns the center point of this Selectable */
	public Point getSelectionPoint();
	/** Returns the selection points of this Selectable */
	public Point[] getSelectionPoints();
	/** Returns the center point of this selectable */
	public Point getCenterPoint();
	
	/** Returns the distance of given point to this Selectable */
	public int getDistance(Point p);
	
	/** Returns true if this Selectable should be selectable */
	public boolean isSelectable();
}