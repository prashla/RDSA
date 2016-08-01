
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
import gld.sim.*;
import gld.infra.*;

import java.awt.*;
import java.awt.event.*;

/**
 * Tool to set spawning and destination frequencies.
 *
 * A little HOW-TO:
 * - Open infrastructure in GLDsim.
 * - Select 'Edge config' tool.
 * - Click on an EdgeNode to select it.
 * - Right-click anywhere on the view to deselect.
 * - Select a roaduser type from the dropdown box in the toolbar if you like.
 * - To change spawning frequency (for currently selected roaduser type):
 *		- Left-click on currently selected EdgeNode (with the magenta rectangle around it).
 *		- Hold the mousebutton and drag left to decrease, right to increase the spawning frequency.
 * - To change destination frequencies (for currently selected roaduser type):
 *		- When an EdgeNode is selected, you can change the chance for each other EdgeNode that
 *			a newly spawned Roaduser will go there.
 *		- Left-click on the EdgeNode you want to change the frequency for.
 *		- Hold the mousebutton and drag left to decrease, right to increase the destination frequency.
 *
 * // This really works a *lot* easier than any dialog I could think of. :) //
 * 
 * @author Group GUI
 * @version 1.0
 */

public class EdgeNodeTool implements Tool
{
	RoaduserPanel roaduserPanel;
	Choice roaduserChoice;
	
	EdgeNode selected = null, modify = null;
	Point lastPoint;
	
	EdgeNode[] edges;

	Controller controller;
	View view;
	
	int curRuType;

	/** 
	* Creates a new EdgeNodeTool.
	* @param c The SimController this tool is used in.
	* @param _view The View this tool is used on.
	* @param _edges All EdgeNodes in the current infrastructure.
	*/
	public EdgeNodeTool(Controller c) {
		controller = c;
		view = controller.getView();
		edges = controller.getModel().getInfrastructure().getEdgeNodes_();
		roaduserPanel = new RoaduserPanel();
	}
	
	/**
	* Sets the current set of EdgeNodes.
	*/
	public void setEdgeNodes(EdgeNode[] _edges)
	{
		edges = _edges;
		selected = null;
		modify = null;
	}

	/**
	* Sets the current roaduser type.
	* @param sel Description of the roaduser type to set as current.
	*/
	protected void setRuType(String sel)
	{
		curRuType = RoaduserFactory.getTypeByDesc(sel);
		selected = null;
		modify = null;
		view.repaint();
	}

	public void mousePressed(View view, Point p, Tool.Mask mask)
	{
		if (selected != null) { // EdgeNode selected
			if (mask.isRight()) { // right-click
				selected = null; // deselect
				modify = null;
			}
			else if(mask.isLeft()) // left-click
			{
				lastPoint = p;
				modify = findEdgeNode(p); // modify frequency for 'modify'
			}
			view.repaint();
			return;
		}
		if(mask.isRight()) {
			Node node = findNode(p);
			if(node!=null)   
				((SimController)controller).getSimModel().getTLController().trackNode(node.getId()); 

		}
		if (!mask.isLeft()) return;
		
		selected = findEdgeNode(p); // select a new EdgeNode
		view.repaint();			
	}
	
	/** Returns the EdgeNode at given Point. */
	protected EdgeNode findEdgeNode(Point p)
	{
		EdgeNode best = null;
		int dist = Integer.MAX_VALUE;
		
		// Find edgenode that was clicked.
		for (int i=0; i < edges.length; i++) {
			if (edges[i].getBounds().contains(p) && edges[i].getDistance(p) < dist) {
				best = edges[i];
				dist = edges[i].getDistance(p);
			}
		}
		return best;
	}

	/** Returns the EdgeNode at given Point. */
	protected Node findNode(Point p)
	{
		Node best = null;
		int dist = Integer.MAX_VALUE;
		
		Node[] allNodes = controller.getModel().getInfrastructure().getAllNodes();
		int num_nodes = allNodes.length;
		
		// Find edgenode that was clicked.
		for (int i=0; i < num_nodes; i++) {
			if (allNodes[i].getBounds().contains(p) && allNodes[i].getDistance(p) < dist) {
				best = allNodes[i];
				dist = allNodes[i].getDistance(p);
			}
		}
		return best;
	}

	
	public void mouseReleased(View view, Point p, Tool.Mask mask)
	{
		if (modify != null && modify == selected) // currently dragging
		{
			SimModel sm = (SimModel)controller.getModel();
			float freqChange = (p.x - lastPoint.x) / (float)250.0;
			float spawn = modify.getSpawnFrequency(curRuType);
			spawn += freqChange;
			if(spawn > 1) spawn = (float)1.0;
			if(spawn < 0) spawn = (float)0.0;
			sm.setSpawnFrequency(modify, curRuType, spawn);
		}
		modify = null; // no longer modifying a frequency
		return;
	}

	public void mouseMoved(View view, Point p, Tool.Mask mask) 
	{
		if(modify != null) // currently dragging
		{
			SimModel sm = (SimModel)controller.getModel();
			float freqChange = (p.x - lastPoint.x) / (float)250.0;
			if(modify == selected) // modifying spawn frequency
			{
				float spawn = modify.getSpawnFrequency(curRuType);
				spawn += freqChange;
				if(spawn > 1) spawn = (float)1.0;
				if(spawn < 0) spawn = (float)0.0;
				modify.setSpawnFrequency(curRuType, spawn);
			}
			else // modifying destination frequency
			{
				float dest = selected.getDestFrequency(modify.getId(), curRuType);
				dest += freqChange;
				if(dest > 1) dest = (float)1.0;
				if(dest < 0) dest = (float)0.0;
				selected.setDestFrequency(modify.getId(), curRuType, dest);
			}
			lastPoint = p;
			view.repaint();
		}
	}
	
	public int overlayType() { return 1; }

	public void paint(Graphics g) throws GLDException {
		if(selected == null) return;
		Rectangle r = selected.getBounds();
		r.grow(3, 3);		
		g.setPaintMode();
		g.setColor(Color.magenta);
		g.drawRect(r.x, r.y, r.width, r.height);

		EdgeNode curEdge;
		for(int i=0; i<edges.length; i++)
		{
			curEdge = edges[i];
			r = curEdge.getBounds();
			
			if(curEdge != selected)
			{
				String dest = "" + selected.getDestFrequency(curEdge.getId(), curRuType);
				if(dest.length() > 5) dest = dest.substring(0,5);
				g.drawString(dest, r.x, r.y - 8);
			}	
			else 
			{
				String spawn = "" + selected.getSpawnFrequency(curRuType);
				if(spawn.length() > 5) spawn = spawn.substring(0,5);
				g.drawString(spawn, r.x, r.y - 8);
			}
		}
	}
	
	public Panel getPanel() { return roaduserPanel; }
	
	/** Panel containing a dropdown box for all concrete Roaduser types. */
	protected class RoaduserPanel extends Panel
	{
		public RoaduserPanel()
		{
			setLayout(null);
			
			roaduserChoice = new Choice();	
			String[] descs;
			if(selected!=null) {
				SpawnFrequency[] freqs = modify.getSpawnFrequencies();
				descs = new String[freqs.length];
				for(int i=0;i<descs.length;i++)
					descs[i] = RoaduserFactory.getDescByType(freqs[i].ruType);
			}
			else
			{
				descs = RoaduserFactory.getConcreteTypeDescs();
			}
			
			for(int i=0; i<descs.length; i++)
				roaduserChoice.add(descs[i]);
			roaduserChoice.addItemListener(new MyListener());
			this.add(roaduserChoice);
			setRuType(descs[0]);
			
			roaduserChoice.setBounds(0, 0, 100, 24);
			setSize(100, 24);
		}
	
		protected class MyListener implements ItemListener
		{
			public void itemStateChanged(ItemEvent e)
			{
				String sel = ((Choice)e.getSource()).getSelectedItem();
				setRuType(sel);
			}
		}	
	}
}