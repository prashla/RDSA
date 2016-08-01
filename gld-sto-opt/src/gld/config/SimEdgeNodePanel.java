
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
import java.util.*;

import gld.*;
import gld.sim.*;
import gld.infra.*;
import gld.utils.*;
import gld.sim.stats.*;

/**
 *
 * @author Group GUI
 * @version 1.0
 */

public class SimEdgeNodePanel extends ConfigPanel implements ItemListener, ActionListener
{
	EdgeNode edgenode;

	TextField spawnFreq;
	Choice spawnTypes;
	Button setSpawn;
	
	Hyperlink wqlLink, twtLink, ruaLink, roadLink, nodeLink;
	Label[] queue;
	
	public SimEdgeNodePanel(ConfigDialog cd, EdgeNode e) {
		super(cd);

		Label clab = new Label("Connects:");
		clab.setBounds(0, 0, 100, 20);
		add(clab);
		
		roadLink = new Hyperlink();
		roadLink.addActionListener(this);
		roadLink.setBounds(100, 0, 100, 20);
		add(roadLink);
		
		Label wlab = new Label("With:");
		wlab.setBounds(0, 20, 100, 20);
		add(wlab);

		nodeLink = new Hyperlink();
		nodeLink.addActionListener(this);
		nodeLink.setBounds(100, 20, 100, 20);
		add(nodeLink);


		wqlLink = new Hyperlink("Track waiting queue length");
		wqlLink.addActionListener(this);
		wqlLink.setBounds(0, 50, 200, 20);
		add(wqlLink);

		twtLink = new Hyperlink("Track trip waiting time");
		twtLink.addActionListener(this);
		twtLink.setBounds(0, 70, 200, 20);
		add(twtLink);
		
		ruaLink = new Hyperlink("Track roadusers arrived");
		ruaLink.addActionListener(this);
		ruaLink.setBounds(0, 90, 200, 20);
		add(ruaLink);




		Label lab = new Label("Spawnfrequency for");
		lab.setBounds(0, 120, 120, 20);
		add(lab);

		spawnTypes = new Choice();
		spawnTypes.addItemListener(this);

		String[] descs = RoaduserFactory.getConcreteTypeDescs();
		for (int i=0; i < descs.length; i++)
			spawnTypes.addItem(descs[i]);

		spawnTypes.setBounds(0, 140, 100, 20);
		add(spawnTypes);
		
		lab = new Label("is");
		lab.setBounds(105, 140, 15, 20);
		add(lab);

		spawnFreq = new TextField();
		spawnFreq.setBounds(120, 140, 40, 20);
		spawnFreq.addActionListener(this);
		add(spawnFreq);
		
		setSpawn = new Button("Set");
		setSpawn.addActionListener(this);
		setSpawn.setBounds(170, 140, 50, 20);
		add(setSpawn);
		
		lab = new Label("Waiting in queue:");
		lab.setBounds(200, 0, 150, 20);
		add(lab);

		int nrtypes = RoaduserFactory.statArrayLength();
		queue = new Label[nrtypes];
		for (int i=0; i < nrtypes; i++) {
			lab = new Label();
			lab.setBounds(200, i * 20 + 20, 150, 20);
			add(lab);
			queue[i] = lab;
		}



		setEdgeNode(e);
	}

	public void reset() {
		setSpawnFreq();
		try {
			int nrtypes = RoaduserFactory.statArrayLength();
			int[] nrwaiting = new int[nrtypes];
			Roaduser ru;
			ListIterator li = edgenode.getWaitingQueue().listIterator();
			while (li.hasNext()) {
				ru = (Roaduser)li.next();
				nrwaiting[RoaduserFactory.getStatIndexByType(ru.getType())]++;
				nrwaiting[0]++;
			}
			for (int i=0; i < nrtypes; i++) {
				queue[i].setText(nrwaiting[i] + " - " + RoaduserFactory.getDescByStatIndex(i));
			}
		}
		// SimModel thread changed the queue while we were updating, try again.
		catch (ConcurrentModificationException e) {
			reset();
		}
	}

	public void setSpawnFreq() {
		int type = getSpawnType();
		float freq = edgenode.getSpawnFrequency(type);
		spawnFreq.setText("" + (freq > 0 ? freq : 0));
	}

	public void setSpawnType() {
		try {
			SimModel sm = (SimModel)confd.getController().getModel();
			float fr = Float.parseFloat(spawnFreq.getText());
			sm.setSpawnFrequency(edgenode, getSpawnType(), fr);
		}
		catch (NumberFormatException ex) {
			confd.showError("You must enter a float");
		}
	}
		
	
	public void setEdgeNode(EdgeNode e) {
		edgenode = e;
		confd.setTitle(edgenode.getName());
		reset();
		setSpawnType();
		
		Road road = edgenode.getRoad();
		if (road != null) {
			roadLink.setText(road.getName());
			roadLink.setEnabled(true);
			nodeLink.setText(road.getOtherNode(edgenode).getName());
			nodeLink.setEnabled(true);
		}
		else {
			roadLink.setText("null");
			roadLink.setEnabled(false);
			nodeLink.setText("null");
			nodeLink.setEnabled(false);
		}
	}
	
	public void itemStateChanged(ItemEvent e) {
		Object source = e.getItemSelectable();
		if (source == spawnTypes) setSpawnFreq();
	}

	/** Returns the currently selected roaduser type */	
	public int getSpawnType() {
		int[] types = RoaduserFactory.getConcreteTypes();
		return types[spawnTypes.getSelectedIndex()];
	}

	
	public void actionPerformed(ActionEvent e)
	{
		Object source = e.getSource();
		
		if (source == setSpawn || source == spawnFreq) setSpawnType();
		else if (source == wqlLink) track(TrackerFactory.SPECIAL_QUEUE);
		else if (source == twtLink) track(TrackerFactory.SPECIAL_WAIT);
		else if (source == ruaLink) track(TrackerFactory.SPECIAL_ROADUSERS);
		else if (source == roadLink) confd.selectObject(edgenode.getRoad());
		else if (source == nodeLink) confd.selectObject(edgenode.getRoad().getOtherNode(edgenode));
	}

	public void track(int type) {
		SimController sc = (SimController)confd.getController();
		try {
			TrackerFactory.showTracker(sc.getSimModel(), sc, edgenode, type);
		}
		catch (GLDException ex) {
			Controller.reportError(ex);
		}
	}
}
	
