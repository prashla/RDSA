
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
import java.awt.event.*;

/**
 * This tool allows you to do everything at once.
 * It implements the RoadTool, NodeTool, MoveTool, ScrollTool, ZoomTool and, later, the EdgeNodeTool
 *
 * @author Group GUI
 * @version 1.0
 */


public class TotalEditTool implements Tool
{
	protected EditController controller;
	protected NodeTypeChoice typePanel;
	
	public TotalEditTool(EditController ec) {
		controller = ec;
		typePanel = new NodeTypeChoice();
	}
	
	public void mousePressed(View view, Point p, Tool.Mask mask) { }
	public void mouseReleased(View view, Point p, Tool.Mask mask) { }
	public void mouseMoved(View view, Point p, Tool.Mask mask) { }

	public int overlayType() { return 0; }
	public void paint(Graphics g) throws GLDException { }
	
	public Panel getPanel() { return typePanel; }


  protected class NodeTypeChoice extends Panel implements ItemListener
  {
  	int nodeType = 2;
  
  	public NodeTypeChoice()
  	{
  		super();
	  	Choice nodeTypeSel = new Choice();
  		nodeTypeSel.add("Edge node");
  		nodeTypeSel.add("Traffic lights");
  		nodeTypeSel.add("No signs");
  		nodeTypeSel.select(1);
  		nodeTypeSel.addItemListener(this);
  		this.add(nodeTypeSel);
  	}
  	
  	public int getNodeType() { return nodeType; }
  	public void setNodeType(int type) { nodeType = type; }

		public void itemStateChanged(ItemEvent e) {
			setNodeType(((Choice) e.getSource()).getSelectedIndex() + 1);
		}
  }
}