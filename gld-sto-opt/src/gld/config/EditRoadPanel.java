
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
import gld.infra.*;
import gld.edit.*;
import gld.utils.*;

/**
 *
 * @author Group GUI
 * @version 1.0
 */

public class EditRoadPanel extends ConfigPanel implements ActionListener, ItemListener
{
	Road road;
	Drivelane lane;
	
	Hyperlink alphaLink, betaLink, laneLink;
	
	List alphaList;
	List betaList;
	Button addAlpha;
	Button addBeta;
	
	Choice typeList;
	Checkbox left, ahead, right;
	Button delete;



	public EditRoadPanel(ConfigDialog cd, Road r)
	{
		super(cd);
		
		Label alab = new Label("Lanes to");
		alab.setBounds(0, 0, 60, 20);
		add(alab);
		
		alphaLink = new Hyperlink();
		alphaLink.setBounds(60, 0, 80, 20);
		alphaLink.addActionListener(this);
		add(alphaLink);
		
		Label blab = new Label("Lanes to");
		blab.setBounds(200, 0, 60, 20);
		add(blab);
		
		betaLink = new Hyperlink();
		betaLink.setBounds(260, 0, 80, 20);
		betaLink.addActionListener(this);
		add(betaLink);
		
		alphaList = new List();
		alphaList.setBounds(0, 25, 150, 40);
		alphaList.addItemListener(this);
		add(alphaList);
		
		betaList = new List();
		betaList.setBounds(200, 25, 150, 40);
		betaList.addItemListener(this);
		add(betaList);
		
		addAlpha = new Button();
		addAlpha.addActionListener(this);
		addAlpha.setBounds(0, 70, 150, 25);
		add(addAlpha);
		
		addBeta = new Button();
		addBeta.addActionListener(this);
		addBeta.setBounds(200, 70, 150, 25);
		add(addBeta);




		laneLink = new Hyperlink();
		laneLink.setBounds(0, 130, 100, 20);
		laneLink.addActionListener(this);
		add(laneLink);

		Label allows = new Label("allows:");
		allows.setBounds(100, 130, 50, 20);
		add(allows);
		
		typeList = new Choice();

		String[] descs = RoaduserFactory.getTypeDescs();
		int[] types = RoaduserFactory.getTypes();

		for (int i=0; i < descs.length; i++) typeList.addItem(descs[i]);
		
		typeList.addItemListener(this);
		typeList.setBounds(150, 130, 100, 20);
		add(typeList);

		Label turnlab = new Label("Roadusers are allowed to:");
		turnlab.setBounds(0, 160, 300, 20);
		add(turnlab);
		
		left = new Checkbox("Turn left");
		left.setBounds(0, 185, 100, 20);
		left.addItemListener(this);
		add(left);
		
		ahead = new Checkbox("Go straight ahead");
		ahead.setBounds(100, 185, 150, 20);
		ahead.addItemListener(this);
		add(ahead);
		
		right = new Checkbox("Turn right");
		right.setBounds(250, 185, 100, 20);
		right.addItemListener(this);
		add(right);
		
		delete = new Button("Delete drivelane 0");
		delete.setBounds(0, 215, 150, 24);
		delete.addActionListener(this);
		add(delete);

		setRoad(r);
	}
	
	public void paint(Graphics g) {
		super.paint(g);
		g.setColor(Color.black);
		g.drawLine(0, 110, ConfigDialog.PANEL_WIDTH, 110);
	}
	
	public void deleteLane() {
		if (road.getNumAllLanes() > 1) {
			EditModel em = (EditModel)confd.getController().getModel();
			try { em.remLane(lane); }
			catch (InfraException ex) { Controller.reportError(ex); }
		}
	}
	
	public void addLane(Node to) {
		EditModel em = (EditModel)confd.getController().getModel();
		try {
			em.addLane(new Drivelane(road), road, to);
		}
		catch (InfraException e) {
			Controller.reportError(e);
		}
	}
	
	public void reset() {
		Drivelane[] lanes = road.getAlphaLanes();
		
		if (lanes.length != alphaList.getItemCount()) {
			alphaList.removeAll();
		
			for (int i=0; i < lanes.length; i++)
				alphaList.add(lanes[i].getName());
			if (lanes.length >= 4) addAlpha.setEnabled(false);
			else addAlpha.setEnabled(true);

			if (lanes.length > 0) {
				alphaList.select(0);
				setLane(lanes[0]);
			}
		}

		lanes = road.getBetaLanes();
		
		if (lanes.length != betaList.getItemCount()) {
			betaList.removeAll();

			for (int i=0; i < lanes.length; i++)
				betaList.add(lanes[i].getName());
			if (lanes.length >= 4) addBeta.setEnabled(false);
			else addBeta.setEnabled(true);

			if (lanes.length > 0 && road.getNumAlphaLanes() <= 0) {
				betaList.select(0);
				setLane(lanes[0]);
			}
		}

		if (road.getNumAllLanes() > 1) delete.setEnabled(true);
		else delete.setEnabled(false);
	}
		
	public void setRoad(Road r) {
		road = r;
		confd.setTitle(road.getName());
		
		alphaList.removeAll();
		betaList.removeAll();
		
		Drivelane[] lanes = road.getAlphaLanes();
		
		for (int i=0; i < lanes.length; i++)
			alphaList.add(lanes[i].getName());
		if (lanes.length >= 4) addAlpha.setEnabled(false);
		else addAlpha.setEnabled(true);

		if (lanes.length > 0) {
			alphaList.select(0);
			setLane(lanes[0]);
		}

		lanes = road.getBetaLanes();

		for (int i=0; i < lanes.length; i++)
			betaList.add(lanes[i].getName());
		if (lanes.length >= 4) addBeta.setEnabled(false);
		else addBeta.setEnabled(true);

		if (lanes.length > 0 && road.getNumAlphaLanes() <= 0) {
			betaList.select(0);
			setLane(lanes[0]);
		}

		if (road.getNumAllLanes() > 1) delete.setEnabled(true);
		else delete.setEnabled(false);

		alphaLink.setText(road.getAlphaNode().getName());
		betaLink.setText(road.getBetaNode().getName());
		addAlpha.setLabel("Add lane to " + road.getAlphaNode().getName());
		addBeta.setLabel("Add lane to " + road.getBetaNode().getName());
	}
	
	public void setLane(Drivelane l) {
		lane = l;

		try {
			left.setState(lane.getTarget(0));
			ahead.setState(lane.getTarget(1));
			right.setState(lane.getTarget(2));
		}
		catch (InfraException e) {
			Controller.reportError(e);
		}
		
		typeList.select(RoaduserFactory.getDescByType(lane.getType()));
		laneLink.setText(lane.getName());
		delete.setLabel("Delete " + lane.getName());
	}

	public void actionPerformed(ActionEvent e) {
		Object source = e.getSource();
		if (source == alphaLink)
			confd.selectObject(road.getAlphaNode());
		else if (source == betaLink)
			confd.selectObject(road.getBetaNode());
		else if (source == laneLink)
			confd.selectObject(lane);
		else if (source == delete)
			deleteLane();
		else if (source == addAlpha)
			addLane(road.getAlphaNode());
		else if (source == addBeta)
			addLane(road.getBetaNode());
	}
	
	public void itemStateChanged(ItemEvent e) {
		ItemSelectable es = e.getItemSelectable();
		EditModel em = (EditModel)confd.getController().getModel();
		
		try {
			if (es == alphaList) {
				setLane(road.getAlphaLanes()[alphaList.getSelectedIndex()]);
				betaList.deselect(betaList.getSelectedIndex());
			}
			else if (es == betaList) {
				setLane(road.getBetaLanes()[betaList.getSelectedIndex()]);
				alphaList.deselect(alphaList.getSelectedIndex());
			}
			else if (es == typeList) {
				int[] types = RoaduserFactory.getTypes();
				em.setLaneType(lane, types[typeList.getSelectedIndex()]);
			}
			else if (es == left)
				em.setLaneTarget(lane, 0, left.getState());
			else if (es == ahead)
				em.setLaneTarget(lane, 1, ahead.getState());
			else if (es == right)
				em.setLaneTarget(lane, 2, right.getState());
		}
		catch (GLDException ex) {
			Controller.reportError(ex);
		}
	}
}