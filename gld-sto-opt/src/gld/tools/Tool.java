
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

import gld.Overlay;
import gld.View;

import java.awt.Panel;
import java.awt.Point;
import java.awt.Graphics;
import java.awt.event.InputEvent;

/**
 * All classes implementing this interface can serve as tool to be used in 
 * either simulator, editor or both.
 *
 * @author Group GUI
 * @version 1.0
 */

public interface Tool extends Overlay
{
	/** 
	* Invoked when the user presses a mouse button.
	*
	* @param view The <code>View</code> that the event originates from.
	* @param p The coordinates of the view the mouse cursor was at when the event was generated.
	* @param mask The mask taken directly from the MouseEvent, wrapped in a Mask object
	*/
	public void mousePressed(View view, Point p, Mask mask);
	
	/** 
	* Invoked when the user releases a mouse button.
	*
	* @param view The <code>View</code> that the event originates from.
	* @param p The coordinates of the view the mouse cursor was at when the event was generated.
	* @param mask The mask taken directly from the MouseEvent, wrapped in a Mask object
	*/
	public void mouseReleased(View view, Point p, Mask mask);
	
	/** 
	* Invoked when the user moves the mouse over the <code>View</code>.
	*
	* @param view The <code>View</code> that the event originates from.
	* @param p The coordinates of the view the mouse cursor was at when the event was generated.
	* @param mask The mask taken directly from the MouseEvent, wrapped in a Mask object
	*/	
	public void mouseMoved(View view, Point p, Mask mask);
	
	/**
	* Invoked when this tool is selected as the current tool.
	* The panel will returned will be added to the end of the toolbar, allowing
	* the current tool to add tool specific controls to the toolbar.
	* The tool should take care of the listeners.
	*/
	public Panel getPanel();


	/**
	* Class used to pass on the modifiers.
	*/
	public static class Mask {
		private int value;
		public Mask(int v) { value = v; }
		public int getValue()						{ return value; }
		public String toString()				{ return "Mask, val: " + value; }
		public boolean isLeft()         { return (value & InputEvent.BUTTON1_MASK) > 0; }
		public boolean isMiddle()       { return (value & InputEvent.BUTTON2_MASK) > 0; }
		public boolean isRight()        { return (value & InputEvent.BUTTON3_MASK) > 0; }
		public boolean isAnyAltDown()   { return isAltDown() || isAltGraphDown(); }
		public boolean isAltDown()      { return (value & InputEvent.ALT_MASK) > 0; }
		public boolean isAltGraphDown() { return (value & InputEvent.ALT_GRAPH_MASK) > 0; }
		public boolean isControlDown()  { return (value & InputEvent.CTRL_MASK) > 0; }
		public boolean isShiftDown()    { return (value & InputEvent.SHIFT_MASK) > 0; }
	}
}