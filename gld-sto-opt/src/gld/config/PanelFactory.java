
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

import java.util.*;

import gld.*;
import gld.infra.*;

/**
 *
 * @author Group GUI
 * @version 1.0
 */

public class PanelFactory
{
	public static final int TYPE_SIM = 1;
	public static final int TYPE_EDIT = 2;
	protected int TYPE;
	
	protected ConfigDialog confd;
	
	protected GeneralPanel gp = null;
	protected NetTunnelPanel ntp = null;
	protected EditDrivelanePanel edp = null;
	protected EditEdgeNodePanel eep = null;
	protected EditJunctionPanel ejp = null;
	protected EditRoadPanel erp = null;
	
	protected SimDrivelanePanel sdp = null;
	protected SimEdgeNodePanel sep = null;
	protected SimJunctionPanel sjp = null;
	protected SimRoadPanel srp = null;

	public PanelFactory(ConfigDialog cd, int t) {
		confd = cd;
		TYPE = t;
	}

	public ConfigPanel createPanel(Selection s) throws ConfigException
	{
		if (!s.isEmpty()) {
			Object firstObject = s.getSelectedObjects().getFirst();
			if (firstObject instanceof Node) return getNodePanel((Node)firstObject);
			if (firstObject instanceof Road) return getRoadPanel((Road)firstObject);
			if (firstObject instanceof Drivelane) return getDrivelanePanel((Drivelane)firstObject);
			throw new ConfigException("Unknown object type");
		}
		if (gp == null) {
			gp = new GeneralPanel(confd);
		}
		gp.reset();
		return gp;
	}





	private ConfigPanel getNodePanel(Node n) throws ConfigException
	{
		
		if (n instanceof EdgeNode) return getEdgeNodePanel((EdgeNode)n);
		if (n instanceof Junction) return getJunctionPanel((Junction)n);
		if (n instanceof NetTunnel) return getNetTunnelPanel((NetTunnel)n);
		throw new ConfigException("Unknown object type");
	}





	private ConfigPanel getRoadPanel(Road r) throws ConfigException
	{
		if (TYPE == TYPE_SIM) return createSRP(r);
		if (TYPE == TYPE_EDIT) return createERP(r);
		throw new ConfigException("Unknown object type");
	}
	
	private ConfigPanel getDrivelanePanel(Drivelane l) throws ConfigException
	{
		if (TYPE == TYPE_SIM) return createSDP(l);
		if (TYPE == TYPE_EDIT) return createEDP(l);
		throw new ConfigException("Unknown object type");
	}
	
	private ConfigPanel getEdgeNodePanel(EdgeNode e) throws ConfigException
	{
		if (TYPE == TYPE_SIM) return createSEP(e);
		if (TYPE == TYPE_EDIT) return createEEP(e);
		throw new ConfigException("Unknown object type");
	}
	
	private ConfigPanel getJunctionPanel(Junction j) throws ConfigException
	{
		if (TYPE == TYPE_SIM) return createSJP(j);
		if (TYPE == TYPE_EDIT) return createEJP(j);
		throw new ConfigException("Unknown object type");
	}
	
	private ConfigPanel getNetTunnelPanel(NetTunnel n) throws ConfigException
	{	return createNTP(n);
	}
	














	private ConfigPanel createSRP(Road r) { 
		if (srp == null) srp = new SimRoadPanel(confd, r);
		else srp.setRoad(r);
		return srp;
	}
	private ConfigPanel createERP(Road r) {
		if (erp == null) erp = new EditRoadPanel(confd, r);
		else erp.setRoad(r);
		return erp;
	}

	private ConfigPanel createSDP(Drivelane l) {
		if (sdp == null) sdp = new SimDrivelanePanel(confd, l);
		else sdp.setLane(l);
		return sdp;
	}
	private ConfigPanel createEDP(Drivelane l) {
		if (edp == null) edp = new EditDrivelanePanel(confd, l);
		else edp.setLane(l);
		return edp;
	}

	private ConfigPanel createSEP(EdgeNode e) {
		if (sep == null) sep = new SimEdgeNodePanel(confd, e);
		else sep.setEdgeNode(e);
		return sep;
	}
	private ConfigPanel createEEP(EdgeNode e) {
		if (eep == null) eep = new EditEdgeNodePanel(confd, e);
		else eep.setEdgeNode(e);
		return eep;
	}

	private ConfigPanel createSJP(Junction j) {
		if (sjp == null) sjp = new SimJunctionPanel(confd, j);
		else sjp.setJunction(j);
		return sjp;
	}
	
	private ConfigPanel createEJP(Junction j) {
		if (ejp == null) ejp = new EditJunctionPanel(confd, j);
		else ejp.setJunction(j);
		return ejp;
	}
	
	private ConfigPanel createNTP(NetTunnel n)
	{	if (ntp == null) ntp = new NetTunnelPanel(confd, n);
		else ntp.setNetTunnel(n);
		return ntp;
	}
}