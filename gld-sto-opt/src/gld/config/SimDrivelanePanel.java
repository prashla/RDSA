
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
import java.util.ListIterator;
import java.util.ConcurrentModificationException;
import java.util.Vector;
import java.awt.event.*;

import gld.*;
import gld.infra.*;
import gld.utils.*;
import gld.edit.*;


/**
 *
 * @author Group GUI
 * @version 1.0
 */

public class SimDrivelanePanel extends ConfigPanel implements ActionListener, ItemListener
{
	Drivelane lane;

	Hyperlink alphaLink, betaLink, roadLink;

	Label sign, allows;

	Queue queue;
	Checkbox queueType;
	ScrollPane sp;


	public SimDrivelanePanel(ConfigDialog cd, Drivelane l) {
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

		allows = new Label();
		allows.setBounds(0, 90, 200, 20);
		add(allows);
		
		queue = new Queue();
		sp = new ScrollPane(ScrollPane.SCROLLBARS_ALWAYS);
		sp.setBounds(210, 0, 150, 200);
		sp.add(queue);
		add(sp);

		queueType = new Checkbox("Show free spaces");
		queueType.addItemListener(this);
		queueType.setBounds(210, 210, 150, 20);
		add(queueType);



		setLane(l);
	}
	
	private class Queue extends Panel
	{
		Label[] labels = { };
		Hyperlink[] links = { };
		int counter; // internal counter used by addString and addLink
		public Queue() {
			setLayout(null);
		}
		
		public synchronized void setLane(Drivelane lane) {
			removeAll();
			
			int len = lane.getCompleteLength();
			links = new Hyperlink[len];
			labels = new Label[len];
			Hyperlink link;
			Label lab;
			
			for (int i=0; i < len; i++) {
				link = new Hyperlink();
				link.addActionListener(new RULinkListener(null));
				link.setEnabled(false);
				link.setVisible(false);
				link.setBounds(0, i * 20, 130, 20);
				add(link); 
				links[i] = link;
				
				lab = new Label();
				lab.setBounds(0, i * 20, 130, 20);
				lab.setVisible(false);
				add(lab);
				labels[i] = lab;
			}
			reset();
		}

		public synchronized void reset() {
			int resetTryCount = 10;
			boolean done = false;
			
			while (resetTryCount>0 && !done) try
			{
				boolean sf = queueType.getState();
				ListIterator li = lane.getQueue().listIterator();
				int pos = 0;
				Roaduser ru;
				counter = 0;

				while (li.hasNext())
				{
					ru = (Roaduser)li.next();
					if (sf) {
						while (ru.getPosition() > pos) {
							addString(pos + ": Free block");
							pos++;
						} // while
					} // if
					addLink(ru.getPosition() + ": " + ru.getName(), ru);
					pos += ru.getLength();
					if (sf) {
						for (int i=1; i < ru.getLength(); i++)
							addString(pos + ": -");
					} // if
				} // while
				
				if (pos < lane.getCompleteLength() && sf) {
					while (pos < lane.getCompleteLength()) {
						addString(pos + ": Free block");
						pos++;
					} // while
				} // if

				for (int i=counter; i < labels.length; i++) {
					labels[i].setVisible(false);
					labels[i].repaint();
					links[i].setEnabled(false);
					links[i].setVisible(false);
					links[i].repaint();
				} // for
				done = true;
			} // try
			
			// SimModel thread changed the queue while we were updating, try again.
			catch (ConcurrentModificationException e) {
				done = false;
				resetTryCount--;
				reset();
			}
			catch (NullPointerException e) {
				done = false;
				resetTryCount--;
				reset();
			}
			doLayout();
		}
		
		private void addString(String text) {
			Label lab = labels[counter];
			lab.setVisible(true);
			links[counter].setVisible(false);
			links[counter].setEnabled(false);
			if (!lab.getText().equals(text)) lab.setText(text);
			counter++;
		}
		
		private void addLink(String text, Roaduser ru) {
			Hyperlink link = links[counter];
			link.setVisible(true);
			link.setEnabled(true);
			labels[counter].setVisible(false);

			if (!link.getText().equals(text)) link.setText(text);
			((RULinkListener)link.getActionListeners()[0]).setRoaduser(ru);
			
			counter++;
		}

		private class RULinkListener implements ActionListener
		{
			Roaduser ru;
			public RULinkListener(Roaduser ru) { this.ru = ru; }
			public void setRoaduser(Roaduser ru) { this.ru = ru; }
			public void actionPerformed(ActionEvent e) { confd.showRoaduser(ru); }
		}
		
		public Dimension getPreferredSize() { return new Dimension(130, counter * 20); }
	}




	public void reset() {
		if (lane.getSign().getType() == Sign.NO_SIGN)
			sign.setText("Drivelane has no trafficlight");
		else
			sign.setText(lane.getNumRoadusersWaiting() + " waiting for trafficlight");


		queue.reset();
		sp.doLayout();
	}

	public void setLane(Drivelane l) {
		lane = l;
		confd.setTitle(lane.getName());
		queue.setLane(lane);
		reset();
		
		alphaLink.setText(lane.getNodeLeadsTo().getName());
		betaLink.setText(lane.getNodeComesFrom().getName());
		roadLink.setText(lane.getRoad().getName());
		
		allows.setText("Drivelane allows " + RoaduserFactory.getDescByType(lane.getType()));
	}

	public void itemStateChanged(ItemEvent e) {
		reset();
	}

	public void actionPerformed(ActionEvent e) {
		Object source = e.getSource();
		if (source == alphaLink)
			confd.selectObject(lane.getNodeLeadsTo());
		else if (source == betaLink)
			confd.selectObject(lane.getNodeComesFrom());
		else if (source == roadLink)
			confd.selectObject(lane.getRoad());
	}
}