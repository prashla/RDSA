
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

package gld.sim;

import gld.*;
import gld.infra.*;
import gld.sim.stats.TrackerFactory;

import java.awt.MenuItem;
import java.awt.MenuShortcut;
import java.awt.PopupMenu;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.lang.Exception;


/**
 *
 * Factory for creating popup menus for editor
 *
 * @author Group GUI
 * @version 1.0
 */

public class SimPopupMenuFactory
{
	protected SimController controller;
	
	public SimPopupMenuFactory(SimController con) {
		controller = con;
	}
	
	/**
	 * Creates a right-click PopupMenu for the given object.
	 * A listener is added to the menu as well.
	 */
	public PopupMenu getPopupMenuFor(Selectable obj) throws PopupException
	{
		if (obj instanceof Node) return getNodeMenu((Node)obj);
		if (obj instanceof Road) return getRoadMenu((Road)obj);
		if (obj instanceof Drivelane) return getDrivelaneMenu((Drivelane)obj);
		throw new PopupException("Unknown object type");
	}


	/* Node popup menu's & listeners */

	protected PopupMenu getNodeMenu(Node n) throws PopupException
	{
		PopupMenuListener pml = null;
		if (n instanceof EdgeNode) return getEdgeNodeMenu((EdgeNode)n);
		if (n instanceof Junction) return getJunctionMenu((Junction)n);
		if (n instanceof NetTunnel) return getNetTunnelMenu((NetTunnel)n);
		throw new PopupException("Unknown Node type");
	}

	// EdgeNode popup menu
	
	protected PopupMenu getEdgeNodeMenu(EdgeNode edge)
	{
		PopupMenu menu = new PopupMenu();
		PopupMenuListener pml = new EdgeNodePopupListener(edge, controller);

		String[] items = { "Track waiting queue length", "Track roadusers arrived", "Track trip waiting time" };
		MenuItem item;

		for(int i=0; i<items.length; i++)
		{
	  	item = new MenuItem(items[i]);
  		item.addActionListener(pml);
  		menu.add(item);
		}

		menu.add(new MenuItem("-"));
		
  	item = new MenuItem("Properties...", new MenuShortcut(KeyEvent.VK_ENTER));
  	item.addActionListener(pml);
  	menu.add(item);

  	return menu;
	}
	

	protected static class EdgeNodePopupListener implements PopupMenuListener
	{
		SimController controller;
		EdgeNode node;
		
		public EdgeNodePopupListener(EdgeNode n, SimController con) {
			controller = con;
			node = n;
		}

		public void actionPerformed(ActionEvent e) {
			String sel = e.getActionCommand();
			try
			{
				if (sel.equals("Properties..."))
					controller.showConfigDialog();
				else if(sel.equals("Track waiting queue length"))
					TrackerFactory.showTracker(controller.getSimModel(), controller, node, TrackerFactory.SPECIAL_QUEUE);
				else if(sel.equals("Track trip waiting time"))
					TrackerFactory.showTracker(controller.getSimModel(), controller, node, TrackerFactory.SPECIAL_WAIT);
				else if(sel.equals("Track roadusers arrived"))
					TrackerFactory.showTracker(controller.getSimModel(), controller, node, TrackerFactory.SPECIAL_ROADUSERS);
			}
			catch(GLDException exc) { controller.showError(exc.toString()); }
		}
	}
	
	// NetTunnel popup menu
	
	protected PopupMenu getNetTunnelMenu(NetTunnel tunnel)
	{
		PopupMenu menu = new PopupMenu();
		PopupMenuListener pml = new NetTunnelPopupListener(tunnel, controller);

		String[] items = { "Track roadusers arrived", "Track trip waiting time",
								 "Track waiting (receive) queue","Track send queue"};
		MenuItem item;

		for(int i=0; i<items.length; i++)
		{
	  	item = new MenuItem(items[i]);
  		item.addActionListener(pml);
  		menu.add(item);
		}

		menu.add(new MenuItem("-"));
		
  	item = new MenuItem("Properties...", new MenuShortcut(KeyEvent.VK_ENTER));
  	item.addActionListener(pml);
  	menu.add(item);

  	return menu;
	}
	
	protected static class NetTunnelPopupListener implements PopupMenuListener
	{	SimController controller;
		NetTunnel node;
		
		public NetTunnelPopupListener(NetTunnel n, SimController con) {
			controller = con;
			node = n;
		}

		public void actionPerformed(ActionEvent e) {
			String sel = e.getActionCommand();
			try
			{
				if (sel.equals("Properties..."))
					controller.showConfigDialog();
				else if(sel.equals("Track trip waiting time"))
					TrackerFactory.showTracker(controller.getSimModel(), controller, node, TrackerFactory.SPECIAL_WAIT);
				else if(sel.equals("Track roadusers arrived"))
					TrackerFactory.showTracker(controller.getSimModel(), controller, node, TrackerFactory.SPECIAL_ROADUSERS);
				else if(sel.equals("Track waiting (receive) queue"))
					TrackerFactory.showTracker(controller.getSimModel(),controller,node, TrackerFactory.SPECIAL_QUEUE);
				else if(sel.equals("Track send queue"))
					TrackerFactory.showTracker(controller.getSimModel(), controller,node, TrackerFactory.NETTUNNEL_SEND);
					
				
			}
			catch(GLDException exc) { controller.showError(exc.toString()); }
		}
	}
	

	// Junction popup menu

	protected PopupMenu getJunctionMenu(Junction junction)
	{
		PopupMenu menu = new PopupMenu();
		PopupMenuListener pml = new JunctionPopupListener(junction, controller);

		String[] items = { "Track roadusers that crossed", "Track junction waiting time" };
		MenuItem item;

		for(int i=0; i<items.length; i++)
		{
	  	item = new MenuItem(items[i]);
  		item.addActionListener(pml);
  		menu.add(item);
		}

		menu.add(new MenuItem("-"));
		
  	item = new MenuItem("Properties...", new MenuShortcut(KeyEvent.VK_ENTER));
  	item.addActionListener(pml);
  	menu.add(item);

  	return menu;
	}

	protected static class JunctionPopupListener implements PopupMenuListener
	{
		SimController controller;
		Junction node;
		
		public JunctionPopupListener(Junction n, SimController con) {
			controller = con;
			node = n;
		}
	
		public void actionPerformed(ActionEvent e) {
			String sel = e.getActionCommand();
			if (sel.equals("Properties..."))
				controller.showConfigDialog();
			else
				try { 
					if(sel.equals("Track roadusers that crossed"))
						TrackerFactory.showTracker(controller.getSimModel(), controller, node, TrackerFactory.JUNCTION_ROADUSERS); 
					else if(sel.equals("Track junction waiting time"))
						TrackerFactory.showTracker(controller.getSimModel(), controller, node, TrackerFactory.JUNCTION_WAIT); 
				}
				catch(GLDException exc) {}
		}
	}



	/* Road popup menu & listeners */

	protected PopupMenu getRoadMenu(Road r)
	{
		PopupMenu menu = new PopupMenu();
		PopupMenuListener pml = new RoadPopupListener(r, controller);

  	MenuItem item = new MenuItem("Properties...", new MenuShortcut(KeyEvent.VK_ENTER));
  	item.addActionListener(pml);
  	menu.add(item);

  	return menu;
	}
	
	protected static class RoadPopupListener implements PopupMenuListener
	{
		SimController controller;
		Road road;
		
		public RoadPopupListener(Road r, SimController con) {
			controller = con;
			road = r;
		}
	
		public void actionPerformed(ActionEvent e) {
			if (e.getActionCommand().equals("Properties..."))
				controller.showConfigDialog();
		}
	}

	
	/* Drivelane popup menu & listeners */
	
	protected PopupMenu getDrivelaneMenu(Drivelane l)
	{
		PopupMenu menu = new PopupMenu();
		PopupMenuListener pml = new LanePopupListener(l, controller);

  	MenuItem item = new MenuItem("Properties...", new MenuShortcut(KeyEvent.VK_ENTER));
  	item.addActionListener(pml);
  	menu.add(item);

  	return menu;
	}
	
	protected static class LanePopupListener implements PopupMenuListener
	{
		SimController controller;
		Drivelane lane;
		
		public LanePopupListener(Drivelane l, SimController con) {
			controller = con;
			lane = l;
		}
	
		public void actionPerformed(ActionEvent e) {
			if (e.getActionCommand().equals("Properties..."))
				controller.showConfigDialog();
		}
	}





	/* Popup menu listener interface */

	protected static interface PopupMenuListener extends ActionListener { }
}