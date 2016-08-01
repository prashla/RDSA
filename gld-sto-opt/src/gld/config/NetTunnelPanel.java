
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
import gld.infra.*;
import gld.sim.SimulationRunningException;
import gld.utils.*;

/**
 *
 * @author Siets El snel & J Moritz
 */
 
public class NetTunnelPanel extends ConfigPanel implements ActionListener
{
	NetTunnel netTunnel;
	
	TextField lpField, rpField, rhField;
	Button lpSet, rpSet, rhSet;
	Hyperlink roadLink, nodeLink;
	
	public NetTunnelPanel(ConfigDialog cd, NetTunnel netTunnel) {
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


		Label lplab = new Label("Local port:");
		lplab.setBounds(0, 50, 100, 20);
		add(lplab);
		
		lpField = new TextField();
		lpField.addActionListener(this);
		lpField.setBounds(100, 50, 100, 20);
		add(lpField);
		
		lpSet = new Button("Set");
		lpSet.addActionListener(this);
		lpSet.setBounds(210, 50, 40, 20);
		add(lpSet);


		Label rplab = new Label("Remote port:");
		rplab.setBounds(0, 75, 100, 20);
		add(rplab);

		rpField = new TextField();
		rpField.addActionListener(this);
		rpField.setBounds(100, 75, 100, 20);
		add(rpField);

		rpSet = new Button("Set");
		rpSet.addActionListener(this);
		rpSet.setBounds(210, 75, 40, 20);
		add(rpSet);


		Label rhlab = new Label("Remote host:");
		rhlab.setBounds(0, 100, 100, 20);
		add(rhlab);
		
		rhField = new TextField();
		rhField.addActionListener(this);
		rhField.setBounds(100, 100, 100, 20);
		add(rhField);

		rhSet = new Button("Set");
		rhSet.addActionListener(this);
		rhSet.setBounds(210, 100, 40, 20);
		add(rhSet);




		setNetTunnel(netTunnel);
	}

	public void reset() {
		Road road = netTunnel.getRoad();
		if (road != null) {
			roadLink.setText(road.getName());
			roadLink.setEnabled(true);
			nodeLink.setText(road.getOtherNode(netTunnel).getName());
			nodeLink.setEnabled(true);
		}
		else {
			roadLink.setText("null");
			roadLink.setEnabled(false);
			nodeLink.setText("null");
			nodeLink.setEnabled(false);
		}
	}

	public void setNetTunnel(NetTunnel nt) {
		netTunnel = nt;
		confd.setTitle(netTunnel.getName());
		reset();
		
		lpField.setText(netTunnel.getLocalPort()+"");
		rpField.setText(netTunnel.getRemotePort()+"");
		rhField.setText(netTunnel.getRemoteHostname()+"");
	}

	public void ok() {
		setLocalPort();
		setRemotePort();
		setRemoteHost();
	}

	
	public void actionPerformed (ActionEvent e)
	{
		Object source = e.getSource();
		if (source == lpField || source == lpSet) setLocalPort();
		else if (source == rpField || source == rpSet) setRemotePort();
		else if (source == rhField || source == rhSet) setRemoteHost();
		else if (source == roadLink) confd.selectObject(netTunnel.getRoad());
		else if (source == nodeLink) confd.selectObject(netTunnel.getRoad().getOtherNode(netTunnel));
	}


	public void setLocalPort()
	{
		try
		{
			int pn = Integer.parseInt(lpField.getText());
			if (pn > 65535 || pn < 0) {
					confd.showError ("I have set the local port number to "+
						netTunnel.getRemotePort()+", but that is not a legal value if "+
						"you use TCP (you almost certainly do). The tunnel probably "+
						"won't work.");
			}
			
			netTunnel.setLocalPort(pn);
		}
		catch (NumberFormatException e) {
			confd.showError ("You can only use integers as port numbers.");
		}
		catch (SimulationRunningException e) {
			confd.showError("You must stop the sim before changing the "+
						"properties of this NetTunnel");
		}
	}
	
	public void setRemotePort()
	{
		try
		{
			int pn = Integer.parseInt(rpField.getText());
			if (pn > 65535 || pn < 0) {
					confd.showError ("I have set the remote port number to "+
						netTunnel.getRemotePort()+", but that is not a legal value if "+
						"you use TCP (you almost certainly do). The tunnel probably "+
						"won't work.");
			}
			
			netTunnel.setRemotePort(pn);
		}
		catch (NumberFormatException e) {
			confd.showError ("You can only use integers as port numbers.");
		}
		catch (SimulationRunningException e) {
			confd.showError("You must stop the sim before changing the "+
						"properties of this NetTunnel");
		}
	}
	
	public void setRemoteHost()
	{
		try {
			netTunnel.setRemoteHostname(rhField.getText());
		}
		catch (SimulationRunningException e) {
			confd.showError("You must stop the sim before changing the "+
						"properties of this NetTunnel");
		}
	}
}