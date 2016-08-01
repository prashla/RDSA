
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

package gld.sim.stats;

import gld.sim.SimModel;
import gld.sim.SimController;
import gld.infra.Node;
import gld.infra.RoaduserFactory;
import gld.utils.CheckMenu;
import gld.utils.Arrayutils;

import java.awt.*;
import java.awt.event.*;
import java.util.*;

/**
*
* The extended controller for the tracking window, it controls a TrackingView.
* Offers functionality to choose between allTime average and last 1000 roadusers
* and lets you choose which roaduser types the data are drawn of.
*
* @author Group GUI
* @version 1.0
*/

public class ExtendedTrackingController extends TrackingController
{
	ExtendedTrackingView extView;	
	CheckMenu modeCM, ruTypeCM;

	/**
	* Creates a <code>ExtendedTrackingController</code>.
	*
	* @param _model The <code>SimModel</code> statistics should be read from.
	* @param _controller The parent <code>SimController</code>.
	* @param _extView The <code>ExtendedTrackingView</code> to be shown.
	*/
	public ExtendedTrackingController(SimModel _model, SimController _controller, ExtendedTrackingView _extView)
	{
		super(_model, _controller, _extView);
		extView = _extView;
		if(!extView.useModes()) modeCM.setEnabled(false);
	}
	
	protected void addToOptionsMenu(Menu menu)
	{
		Menu submenu; MenuItem item;
		
		menu.add(new MenuItem("-"));
		
  	ETCMenuListener ml = new ETCMenuListener();

		String[] modes = {"all roadusers", "last " + Node.STAT_NUM_DATA + " roadusers"};
		modeCM = new CheckMenu("Track average of", modes, false);
		modeCM.addItemListener(ml);
		modeCM.select(0);
		menu.add(modeCM);

		String[] ruTypes = {"All roadusers"};
		ruTypes = (String[])Arrayutils.addArray(ruTypes, RoaduserFactory.getConcreteTypeDescs());
		ruTypeCM = new CheckMenu("Roaduser type", ruTypes, true);
		ruTypeCM.addItemListener(ml);
		ruTypeCM.select(0);
		menu.add(ruTypeCM);
 	}

	/**
	* Listens to the extra menu items.
	*/
	public class ETCMenuListener implements ItemListener
	{
		public void itemStateChanged(ItemEvent e) {
			CheckMenu menu = (CheckMenu)e.getItemSelectable();
			if(menu == modeCM)
				extView.setAllTime(modeCM.getSelectedIndex() == 0);
			else {
				CheckboxMenuItem[] citems = ruTypeCM.getItems();
				for(int i=0; i<citems.length; i++) {
					int statIndex = RoaduserFactory.getStatIndexByDesc(citems[i].getLabel());
					extView.showGraph(statIndex, citems[i].getState());
				}
			}
		}
	}
}