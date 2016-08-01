
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
import gld.utils.*;
import gld.algo.dp.*;
import gld.algo.tlc.*;
import gld.sim.stats.*;

import java.awt.*;
import java.util.*;
import java.awt.event.*;

/**
 *
 * The MenuBar for the editor
 *
 * @author Group GUI
 * @version 1.0
 */

public class SimMenuBar extends MenuBar
{
	SimController controller;
	
	SpeedMenu speedMenu;
	TLCMenu tlcMenu;
	DPMenu dpMenu;
	CheckboxMenuItem viewEnabled;
	CheckboxMenuItem cycleCounterEnabled;
	
	public SimMenuBar(SimController sc, String[] speedTexts) {
		
		String[] trackers = {
			"Total waiting queue length",
			"Total roadusers arrived",
			"Average trip waiting time",
			"Average junction waiting time"
		};
		
		String[] dps = DPFactory.getDescriptions();

		controller = sc;
		
		Menu menu; MenuItem item;
		
		add(new FileMenu(controller, false));
		
		/*  Simulation */

		menu = new Menu("Simulation");
		add(menu);
		SimListener simListener = new SimListener();

		item = new MenuItem("Do one step", new MenuShortcut(KeyEvent.VK_D));
		menu.add(item);
		item.addActionListener(simListener);

		item = new MenuItem("Run", new MenuShortcut(KeyEvent.VK_R));
		menu.add(item);
		item.addActionListener(simListener);

		item = new MenuItem("Pause", new MenuShortcut(KeyEvent.VK_U));
		menu.add(item);
		item.addActionListener(simListener);

		item = new MenuItem("Stop", new MenuShortcut(KeyEvent.VK_P));
		menu.add(item);
		item.addActionListener(simListener);

		item = new MenuItem("Run Series", new MenuShortcut(KeyEvent.VK_S));
		menu.add(item);
		item.addActionListener(simListener);

		/* Speed */
				
		speedMenu = new SpeedMenu(speedTexts);
		menu.add(speedMenu);
		
		/* Statistics */

		menu = new Menu("Statistics");
		add(menu);
		StatsListener statsListener = new StatsListener();

		item = new MenuItem("Show statistics", new MenuShortcut(KeyEvent.VK_T));
		menu.add(item);
		item.addActionListener(statsListener);

		Menu submenu = new Menu("Track");

		for(int i=0; i<trackers.length; i++) {
	    	item = new MenuItem(trackers[i]);
  			submenu.add(item);
				item.addActionListener(statsListener);
		}
		menu.add(submenu);

		menu.add(new MenuItem("-"));

		CheckboxMenuItem citem = new CheckboxMenuItem("Toggle in-view statistics", false);
		menu.add(citem);
		citem.addItemListener(statsListener);



		/* Options */

		menu = new Menu("Options");
		add(menu);
		OptionMenuListener ol = new OptionMenuListener();

		viewEnabled = new CheckboxMenuItem("Toggle view", true);
		menu.add(viewEnabled);
		viewEnabled.setName("view");
		viewEnabled.addItemListener(ol);

		cycleCounterEnabled = new CheckboxMenuItem("Toggle cycle counter", true);
		menu.add(cycleCounterEnabled);
		cycleCounterEnabled.setName("cyclecounter");
		cycleCounterEnabled.addItemListener(ol);

		menu.add(new MenuItem("-"));
		
		tlcMenu = new TLCMenu();
		menu.add(tlcMenu);
		
		dpMenu = new DPMenu(dps);
		menu.add(dpMenu);				
				
		menu.add(new MenuItem("-"));

		item = new MenuItem("Open editor", new MenuShortcut(KeyEvent.VK_E));
		menu.add(item);
		item.addActionListener(ol);
		
		item = new MenuItem("Settings...");
		menu.add(item);
		item.addActionListener(ol);
		
		add(new HelpMenu(controller));
	}
		
	public SpeedMenu getSpeedMenu() { return speedMenu; }
	public TLCMenu getTLCMenu() { return tlcMenu; }
	public DPMenu getDPMenu() { return dpMenu; }
		
	protected void setViewEnabled(boolean b) 
	{
		viewEnabled.setState(b);
	}
	
	protected void setCycleCounterEnabled(boolean b)
	{
		cycleCounterEnabled.setState(b);
	}
			
	public class SpeedMenu extends CheckMenu implements ItemListener
	{
		public SpeedMenu(String[] texts) {
			super("Speed", texts);
			addItemListener(this);
			select(1);
		}
			
		public void itemStateChanged(ItemEvent e) {
			controller.setSpeed(getSelectedIndex());
		}
	}



	public class TLCMenu extends Menu implements ItemListener
	{
		CheckMenu[] submenus;
		int selectedMenu = -1;
		
		public TLCMenu()
		{
			super("Traffic light controller");
					
			String[] tlcCats = TLCFactory.getCategoryDescs();
			String[] tlcDescs = TLCFactory.getTLCDescriptions();
			int[][] allTLCs = TLCFactory.getCategoryTLCs();

			submenus = new CheckMenu[tlcCats.length];
			CheckMenu subsubmenu;
			String[] texts;

			for (int i=0; i < tlcCats.length; i++) {
				texts = new String[allTLCs[i].length];
				for (int j=0; j < texts.length; j++)
					texts[j] = tlcDescs[allTLCs[i][j]];
				subsubmenu = new CheckMenu(tlcCats[i], texts);
				add(subsubmenu);
				subsubmenu.addItemListener(this);
				
				submenus[i] = subsubmenu;
			}
			
			setTLC(0, 0);
		}
		
		public void itemStateChanged(ItemEvent e) {
			CheckMenu cm = (CheckMenu)e.getItemSelectable();
			for (int i=0; i < submenus.length; i++) {
				if (submenus[i] == cm) {
					selectedMenu = i;
					controller.setTLC(i, cm.getSelectedIndex());
				}
				else submenus[i].deselectAll();
			}
		}
		
		public void setTLC(int cat, int tlc) {
			for (int i=0; i < submenus.length; i++)
				submenus[i].deselectAll();
			submenus[cat].select(tlc);
		}
		
		public int getCategory() {
			return selectedMenu;
		}
		public int getTLC() {
			return submenus[selectedMenu].getSelectedIndex();
		}
	}
	
	public class DPMenu extends CheckMenu implements ItemListener
	{
		public DPMenu(String[] dps) {
			super("Driving policy", dps);
			addItemListener(this);
			select(0);
		}
		
		public void itemStateChanged(ItemEvent e) {
			controller.setDrivingPolicy(getSelectedIndex());
		}
	}









	/*============================================*/
	/* Listeners                                  */
	/*============================================*/

	/** Listens to the "Statistics" menu */
	public class StatsListener implements ActionListener, ItemListener
	{
		public void actionPerformed(ActionEvent e)
		{
			String sel = e.getActionCommand();

			if(controller.getSimModel().getInfrastructure().getNumNodes() == 0) {
				controller.showError("Please load an infrastructure or simulation before opening any statistics windows.");
				return;
			}

			if (sel.equals("Show statistics")) controller.showStatistics();
			else if (sel.equals("Total waiting queue length")) controller.showTracker(TrackerFactory.TOTAL_QUEUE);
			else if (sel.equals("Total roadusers arrived")) controller.showTracker(TrackerFactory.TOTAL_ROADUSERS);
			else if (sel.equals("Average trip waiting time")) controller.showTracker(TrackerFactory.TOTAL_WAIT);
			else if (sel.equals("Average junction waiting time")) controller.showTracker(TrackerFactory.TOTAL_JUNCTION);
		}

		public void itemStateChanged(ItemEvent e) 
		{
			CheckboxMenuItem item = (CheckboxMenuItem)e.getItemSelectable();
			if(item.getState()) controller.enableOverlay();
			else controller.disableOverlay();
		}
	}


	/** Listens to the "Options" menu */
	public class OptionMenuListener implements ActionListener, ItemListener
	{
		public void actionPerformed(ActionEvent e) {
			String ac = e.getActionCommand();
			if (ac == "Open editor") controller.openEditor();
			else if (ac == "Settings...") controller.showSettings();
		}

		public void itemStateChanged(ItemEvent e) {
			CheckboxMenuItem item = (CheckboxMenuItem)e.getItemSelectable();
			boolean enable = item.getState();
			if(item.getName().equals("view"))
				controller.setViewEnabled(enable);
			else
				controller.setCycleCounterEnabled(enable);
		}
	}


	/** Listens to the "Simulation" menu */
	public class SimListener implements ActionListener
	{
		public void actionPerformed(ActionEvent e)
		{
			String sel = e.getActionCommand();
			if (sel.equals("Do one step")) controller.doStep();
			else if (sel.equals("Run")) controller.unpause();
			else if (sel.equals("Pause")) controller.pause();
			else if (sel.equals("Stop")) controller.stop();
			else if (sel.equals("Run Series")) controller.runSeries();
		}
	}
}