
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

package gld.edit;

import gld.*;
import gld.infra.*;
import gld.edit.*;
import gld.sim.*;
import gld.utils.*;
import gld.xml.*;

import java.awt.Point;
import java.io.IOException;
import java.util.*;

/**
 *
 * The model for the editor.
 *
 * @author Group Model
 * @version 1.0
 */
 
public class EditModel extends Model 
{
	/** Number dispenser for temporary drivelane ids */
	public NumberDispenser laneNumbers;
	/** Number dispenser for temporary road ids */
	public NumberDispenser roadNumbers;

	/** Creates a new model for the editor */
	public EditModel() {
		super();
		laneNumbers = new NumberDispenser();
		roadNumbers = new NumberDispenser();
	}
	
	/** Creates a new model for the editor with the specified infrastructure */
	public EditModel(Infrastructure i) {
		super(i);
		laneNumbers = new NumberDispenser();
		roadNumbers = new NumberDispenser();
	}
	
	/**
	 * Removes a list objects from the infrastructure
	 *
	 * @param list The list of objects to remove
	 */
	public void remObjects(LinkedList list) throws InfraException {
		ListIterator it = list.listIterator();
		Object o;
		while (it.hasNext()) {
			o = it.next();
			if (o instanceof Drivelane) remLane((Drivelane)o);
		}
		it = list.listIterator();
		while (it.hasNext()) {
			o = it.next();
			if (o instanceof Road) remRoad((Road)o);
		}
		it = list.listIterator();
		while (it.hasNext()) {
			o = it.next();
			if (o instanceof Node) remNode((Node)o);
		}
		setChanged();
		notifyObservers();
	}
			
	
	/**
	 * Adds a new node to the infrastructure
	 *
	 * @param coords The pixel coordinates of the new node
	 * @param type The type of the new node
	 */
	public void addNode(Point coords, int type) throws InfraException
	{
		Node n = null;
		switch (type) {
			case 1: { n = new EdgeNode(coords); break; }
			case 2: { n = new Junction(coords); break; }
			case 3: { n = new NonTLJunction(coords); break; }
			case 4: { n= new NetTunnel(coords); break; }
		}
		if (n == null) throw new InfraException("No node with type " + type + " exists");
		infra.addNode(n);
		setChanged();
		notifyObservers();
	}
	


	/**
	 * Adds a new road to the infrastructure
	 *
	 * @param r The road to add
	 * @param alpha The alpha node to connect the road to
	 * @param posa The connection-position to connect the road at the alpha node at
	 * @param beta The beta node to connect the road to
	 * @param posb The connection-position to connect the road at the beta node at
	 */
	public void addRoad(Road r, Node alpha, int posa, Node beta, int posb) throws InfraException
	{
		if (alpha != r.getAlphaNode() && alpha != r.getBetaNode())
			throw new InfraException("Alpha node parameter is not connected to road");
		if (beta != r.getAlphaNode() && beta != r.getBetaNode())
			throw new InfraException("Beta node parameter is not connected to road");

		alpha.addRoad(r, posa);
		beta.addRoad(r, posb);
		
		r.setId(roadNumbers.get());
		Drivelane[] lanes = r.getAlphaLanes();
		for (int i=0; i < lanes.length; i++) lanes[i].setId(laneNumbers.get());
		lanes = r.getBetaLanes();
		for (int i=0; i < lanes.length; i++) lanes[i].setId(laneNumbers.get());
		
		
		updateSigns(alpha);
		updateSigns(beta);

		setChanged();
		notifyObservers();
	}

	/**
	 * Adds a new road to the infrastructure
	 *
	 * @param road The road to add
	 * @param alpha The alpha node to connect the road to
	 * @param posa The connection-position to connect the road at the alpha node at
	 * @param beta The beta node to connect the road to
	 * @param posb The connection-position to connect the road at the beta node at
	 * @param turnPoints The turn points. Angles are calculated.
	 */
	public void addRoad(Road road, Node alpha, int posa, Node beta, int posb, Vector turnPoints) throws InfraException
	{
		if(getConPos(alpha, posa).x != getConPos(beta, posb).x)
		{
    		turnPoints.add(0, getConPos(alpha, posa));
		    turnPoints.add(getConPos(beta, posb));
			Turn[] turns = new Turn[turnPoints.size()-1];
		
			for (int i=0; i < turns.length; i++) {
				Point before = (Point)turnPoints.get(i);
				Point after = (Point)turnPoints.get(i+1);
				
				turns[i] = createTurn(before, avg(before, after), after);
			}
			road.setTurns(turns);
		}
		addRoad(road, alpha, posa, beta,  posb);

		CurveUtils.setupRoadSizes(road,Infrastructure.blockLength);
	}
	
	/** Calculates infra point where a road would connect to node at conpos */
	private Point getConPos(Node node, int conpos)
	{
		if(conpos==0)
		    return new Point(node.getCoord().x, node.getCoord().y - 5 * node.getWidth() - 20);
		else if(conpos==1)
		    return new Point(node.getCoord().x + 5 * node.getWidth() + 20, node.getCoord().y);
		else if(conpos==2)
	            return new Point(node.getCoord().x, node.getCoord().y + 5 * node.getWidth() + 20);
		return new Point(node.getCoord().x - 5 * node.getWidth() - 20, node.getCoord().y);
	}
	
	/**
	 * Create a turn from 3 points.
	 * First is location of turn before it
	 * Second is location of new turn
	 * Third is location of turn after it
	 * The angle of the new Turn lies within 90 and -90 degrees.
	 */
	private Turn createTurn(Point before, Point pos, Point after)
	{
		double dx = after.x - before.x;
		double dy = after.y - before.y;
		double angle = Math.PI * 0.5;
		if (dx != 0) angle = Math.atan(dy / dx);
		else if (dy > 0) angle = -angle;
		if (dx < 0) angle = Math.PI + angle;

		angle = CurveUtils.normalize(-angle);
		return new Turn(pos, angle);
	}
		
	/** Creates new Point that is average of 2 given Points */
	private Point avg(Point a, Point b) {
		return new Point((a.x + b.x) / 2, (a.y + b.y) / 2);
	}




	/**
	 * Adds a drivelane to a road.
	 *
	 * @param lane The drivelane to add to the road.
	 * @param r The road to add the drivelane to.
	 * @param to Roadusers on the drivelane move towards this node.
	 * @throw InfraException If there are 4 drivelanes on that side of the road already.
	 * @throw InfraException If the given road is not connected to the given node.
	 * @throw InfraException If the lane could not be added for any other reason.
	 */
	public void addLane(Drivelane lane, Road r, Node to) throws InfraException
	{
		if (r.getNumInboundLanes(to) >= 4) throw new InfraException("Cannot add more than 4 lanes at one side to a road");
		r.addLane(lane, to);
		lane.setId(laneNumbers.get());
		
		updateSigns(to);
		Node n = r.getOtherNode(to);
		if (n instanceof Junction) ((Junction)n).calculateWidth();
		if (to instanceof Junction) ((Junction)to).calculateWidth();
		
		setChanged();
		notifyObservers();
	}
	
	/**
	 * Sets the type of a drivelane.
	 *
	 * @param lane The drivelane to set the type of.
	 * @param t The new type.
	 * @throw InfraException If the type could not be set.
	 */
	public void setLaneType(Drivelane lane, int t) throws InfraException {
		lane.setType(t);
		setChanged();
		notifyObservers();
	}
	
	/**
	 * Sets or clears a target for a drivelane.
	 * IE the arrow on the road.
	 * Target 0 is left, 1 is straight ahead, 2 is right.
	 *
	 * @param lane The drivelane to set a target of.
	 * @param target The target to set or clear.
	 * @param state True to set the target, false to clear it.
	 * @throw InfraException If the target could not be set.
	 */
	public void setLaneTarget(Drivelane lane, int target, boolean state) throws InfraException {
		lane.setTarget(target, state);
		setChanged();
		notifyObservers();
	}
	
	/**
	 * Moves a node to a new position on the map.
	 * This method will be disabled when turns are implemented.
	 *
	 * @param node The node to move.
	 * @param point The new coordinates of the node.
	 */
	public void moveNode(Node node, Point point) {
		node.setCoord(point);
		setChanged();
		notifyObservers();
	}
	
	/**
	 * Removes a node from the infrastructure
	 *
	 * @param node The node to remove
	 * @throw InfraException If anything bad happens.
	 */
	public void remNode(Node node) throws InfraException
	{
		Road[] roads = node.getAllRoads();
		for (int i=0; i < roads.length; i++)
			if (roads[i] != null)
				remRoad(roads[i]);
		infra.remNode(node);
	}
	
	/**
	 * Removes a road from the infrastructure
	 *
	 * @param road The road to remove from the infrastructure
	 * @throw InfraException If something bad happens.
	 */	
	public void remRoad(Road road) throws InfraException
	{
		Node alpha = road.getAlphaNode();
		Node beta = road.getBetaNode();
		alpha.remRoad(road);
		beta.remRoad(road);
		
		updateSigns(alpha);
		updateSigns(beta);
	}
	
	/**
	 * Removes a drivelane from the infrastructure
	 *
	 * @param lane The drivelane to remove
	 * @throw InfraException If something bad happens.
	 */
	public void remLane(Drivelane lane) throws InfraException
	{
		Road road = lane.getRoad();
		if (road.getNumAllLanes() <= 1) { remRoad(road); return; }
		road.remLane(lane);
		Node alpha = road.getAlphaNode();
		Node beta = road.getBetaNode();
		updateSigns(lane.getSign().getNode());
		if (alpha instanceof Junction) ((Junction)alpha).calculateWidth();
		if (beta instanceof Junction) ((Junction)beta).calculateWidth();

		setChanged();
		notifyObservers();
	}


	private void updateSigns(Node node) throws InfraException {
		updateSigns(node, node.getDesiredSignType());
	}
	
	/**
	 * Makes sure all inbound lanes on given node have correct signs.
	 */
	private void updateSigns(Node node, int type) throws InfraException {
		Drivelane[] lanes = node.getInboundLanes();
		Sign[] signs = new Sign[lanes.length];
		for (int i=0; i < lanes.length; i++) {
			signs[i] = Sign.getInstance(type);
			signs[i].setNode(node);
			signs[i].setLane(lanes[i]);
			lanes[i].setSign(signs[i]);
		}
		node.setSigns(signs);
 	}

	// Specific XMLSerializable implementation
	
	public void load (XMLElement myElement,XMLLoader loader) throws XMLTreeException,IOException,XMLInvalidInputException
	{ 	super.load(myElement,loader);
		loader.load(this,laneNumbers);
		loader.load(this,roadNumbers);
	}
	
	public XMLElement saveSelf ()
	{ 	XMLElement result=super.saveSelf();
	  	result.addAttribute(new XMLAttribute("saved-by","editor"));
	  	return result;
	}
	
	public void saveChilds (XMLSaver saver) throws XMLTreeException,IOException,XMLCannotSaveException
	{ 	super.saveChilds(saver);
		saver.saveObject(laneNumbers);
		saver.saveObject(roadNumbers);
	}
	
}