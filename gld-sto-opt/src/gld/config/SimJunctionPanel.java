
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

package gld.config;

import java.awt.*;
import java.awt.event.*;

import gld.*;
import gld.sim.*;
import gld.sim.stats.*;
import gld.utils.*;
import gld.infra.*;

/**
 *
 * @author Group GUI
 * @version 1.0
 */

public class SimJunctionPanel extends ConfigPanel implements ActionListener, ItemListener
{
	Junction junction;

	Hyperlink[] roadLinks;
	Hyperlink awtLink, rucLink;
	Label nrsigns, width, awtLabel;
	Choice typeList;
	Checkbox awtType;
	
	public SimJunctionPanel(ConfigDialog cd, Junction j) {
		super(cd);
		
		String[] dirs = { "north", "east", "south", "west" };
		
		roadLinks = new Hyperlink[4];
		
		for (int i=0; i < 4; i++) {
			Label lab = new Label("Road " + dirs[i] + ": ");
			lab.setBounds(0, i * 20, 100, 20);
			add(lab);
			
			roadLinks[i] = new Hyperlink();
			roadLinks[i].addActionListener(this);
			roadLinks[i].setBounds(100, i * 20, 100, 20);
			add(roadLinks[i]);
		}
		
		nrsigns = new Label();
		nrsigns.setBounds(200, 0, 200, 20);
		add(nrsigns);
		
		width = new Label();
		width.setBounds(200, 20, 200, 20);
		add(width);


		awtLink = new Hyperlink("Track average waiting time");
		awtLink.addActionListener(this);
		awtLink.setBounds(0, 100, 200, 20);
		add(awtLink);

		rucLink = new Hyperlink("Track roadusers crossed");
		rucLink.addActionListener(this);
		rucLink.setBounds(0, 120, 200, 20);
		add(rucLink);
		
		Label lab = new Label("Average waiting time for");
		lab.setBounds(0, 150, 150, 20);
		add(lab);
		
		typeList = new Choice();
		int nr = RoaduserFactory.statArrayLength();
		for (int i=0; i < nr; i++)
			typeList.add(RoaduserFactory.getDescByStatIndex(i));
		
		typeList.addItemListener(this);
		typeList.setBounds(150, 150, 100, 20);
		add(typeList);
		
		awtLabel = new Label();
		awtLabel.setBounds(260, 150, 100, 20);
		add(awtLabel);
		
		awtType = new Checkbox("Show of last " + Node.STAT_NUM_DATA + " roadusers");
		awtType.addItemListener(this);
		awtType.setBounds(0, 175, 200, 20);
		add(awtType);


		setJunction(j);
	}

	public void reset() {
		awtLabel.setText("is " + junction.getStatistics()[typeList.getSelectedIndex()].getAvgWaitingTime(!awtType.getState()));
	}

	public void setJunction(Junction j) {
		junction = j;
		confd.setTitle(junction.getName());
		reset();

		Road[] roads = junction.getAllRoads();
		
		for (int i=0; i < 4; i++) {
			if (roads[i] != null) {
				roadLinks[i].setText(roads[i].getName());
				roadLinks[i].setEnabled(true);
			}
			else {
				roadLinks[i].setText("null");
				roadLinks[i].setEnabled(false);
			}
		}

		nrsigns.setText("Junction has " + junction.getNumRealSigns() + " trafficlights");
		width.setText("Junction is " + junction.getWidth() + " units wide");
	}

	public void actionPerformed(ActionEvent e)
	{
		SimController sc = (SimController)confd.getController();

		Object source = e.getSource();
		for (int i=0; i < 4; i++)
			if (source == roadLinks[i]) confd.selectObject(junction.getAllRoads()[i]);
		
		try {
			if (source == awtLink)
				TrackerFactory.showTracker(sc.getSimModel(), sc, junction, TrackerFactory.JUNCTION_WAIT);
			else if (source == rucLink)
				TrackerFactory.showTracker(sc.getSimModel(), sc, junction, TrackerFactory.JUNCTION_ROADUSERS);
		}
		catch (GLDException ex) {
			Controller.reportError(ex);
		}
	}
	
	public void itemStateChanged(ItemEvent e) {
		reset();
	}
}