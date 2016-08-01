
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

public class EditDrivelanePanel extends ConfigPanel implements ItemListener, ActionListener
{
	Drivelane lane;

	Hyperlink alphaLink, betaLink, roadLink;

	Label sign;
	
	Choice typeList;
	Checkbox left, ahead, right;
	Button delete;

	public EditDrivelanePanel(ConfigDialog cd, Drivelane l) {
		super(cd);

		Label rlab = new Label("Part of:");
		rlab.setBounds(0, 0, 100, 20);
		add(rlab);
		
		roadLink = new Hyperlink();
		roadLink.addActionListener(this);
		roadLink.setBounds(100, 0, 100, 20);
		add(roadLink);
		
		Label alab = new Label("Leads to:");
		alab.setBounds(0, 20, 100, 20);
		add(alab);
		
		alphaLink = new Hyperlink();
		alphaLink.addActionListener(this);
		alphaLink.setBounds(100, 20, 100, 20);
		add(alphaLink);

		Label blab = new Label("Comes from:");
		blab.setBounds(0, 40, 100, 20);
		add(blab);
		
		betaLink = new Hyperlink();
		betaLink.addActionListener(this);
		betaLink.setBounds(100, 40, 100, 20);
		add(betaLink);

		sign = new Label();
		sign.setBounds(0, 70, 200, 20);
		add(sign);

		Label allows = new Label("This drivelane allows:");
		allows.setBounds(0, 100, 150, 20);
		add(allows);
		
		typeList = new Choice();

		String[] descs = RoaduserFactory.getTypeDescs();
		int[] types = RoaduserFactory.getTypes();

		for (int i=0; i < descs.length; i++) typeList.addItem(descs[i]);
		
		typeList.addItemListener(this);
		typeList.setBounds(150, 100, 100, 20);
		add(typeList);

		Label turnlab = new Label("Roadusers are allowed to:");
		turnlab.setBounds(0, 130, 300, 20);
		add(turnlab);
		
		left = new Checkbox("Turn left");
		left.setBounds(0, 150, 100, 20);
		left.addItemListener(this);
		add(left);
		
		ahead = new Checkbox("Go straight ahead");
		ahead.setBounds(100, 150, 150, 20);
		ahead.addItemListener(this);
		add(ahead);
		
		right = new Checkbox("Turn right");
		right.setBounds(250, 150, 100, 20);
		right.addItemListener(this);
		add(right);
		
		delete = new Button("Delete this drivelane");
		delete.setBounds(0, 190, 150, 24);
		delete.addActionListener(this);
		add(delete);



		setLane(l);
	}
	
	public void reset() {
		if (lane.getSign().getType() == Sign.NO_SIGN)
			sign.setText("Drivelane has no sign");
		else
			sign.setText("Drivelane has normal trafficlight");
		
		if (lane.getRoad().getNumAllLanes() > 1)
			delete.setEnabled(true);
		else
			delete.setEnabled(false);
	}
		
	
	public void setLane(Drivelane l) {
		lane = l;
		confd.setTitle(lane.getName());
		reset();
		
		alphaLink.setText(lane.getNodeLeadsTo().getName());
		betaLink.setText(lane.getNodeComesFrom().getName());
		roadLink.setText(lane.getRoad().getName());

		try {
			left.setState(lane.getTarget(0));
			ahead.setState(lane.getTarget(1));
			right.setState(lane.getTarget(2));
		}
		catch (InfraException e) {
			Controller.reportError(e);
		}
		
		typeList.select(RoaduserFactory.getDescByType(lane.getType()));
	}

	public void itemStateChanged(ItemEvent e) {
		ItemSelectable es = e.getItemSelectable();
		EditModel em = (EditModel)confd.getController().getModel();
		
		try {
			if (es == typeList) {
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

	public void actionPerformed(ActionEvent e) {
		Object source = e.getSource();
		if (source == alphaLink)
			confd.selectObject(lane.getNodeLeadsTo());
		else if (source == betaLink)
			confd.selectObject(lane.getNodeComesFrom());
		else if (source == roadLink)
			confd.selectObject(lane.getRoad());
		else if (source == delete) {
			if (lane.getRoad().getNumAllLanes() > 1) {
				EditModel em = (EditModel)confd.getController().getModel();
				try {
					confd.selectObject(lane.getRoad());
					em.remLane(lane);
				}
				catch (InfraException ex) { Controller.reportError(ex); }
			}
		}
	}
}