
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

package gld.infra;

import gld.*;
import gld.utils.*;
import gld.xml.*;

import java.awt.*;
import java.awt.geom.*;
import java.io.IOException;
import java.util.*;

/**
 *
 * The basic lane.
 *
 * @author Group Datastructures
 * @version 1.0
 */
 
/* Todo: 
 * Possible speedup: cache 'numRoadUsersWaiting' 
 */

public class Drivelane implements XMLSerializable, TwoStageLoader, Selectable
{
    // last road user waiting when the lane is made green (for CDAC algo) 
    public  int lastruStarttime;
    public boolean lastruSetFlag = false;
    
	/** The Id of this drivelane. */
	protected int Id;
	/** The type of drivelane. This is a combination of Roaduser types that may use this drivelane */
	protected int type;
	/** The Road this Drivelane is part of */
	protected Road road;
	/** The roadusers currently on this Drivelane */
	protected LinkedList queue;
	/** The Sign at the end of this Drivelane */
	protected Sign sign;
	/** The directions Roadusers switch lanes to: left, straight ahead and right. */
	protected boolean[] targets = { false, false, false };
	/** The last cycle this Drivelane was asked if it had moved its Roadusers yet */
	protected int cycleAsked;
	/** The last cycle this Drivelane moved its Roadusers */
	protected int cycleMoved;
	/** Data for loading the second stage */
	protected TwoStageLoaderData loadData=new TwoStageLoaderData();
	protected String parentName="model.infrastructure";
	/** A Shape array holding this drivelane's boundaries */
	protected Shape[] bounds = null;

	public Drivelane(Road _road) {
	    road = _road;
	    sign = null;
	    cycleAsked = -1;
	    cycleMoved = -1;
	    Id = -1;
	    type = RoaduserFactory.getTypeByDesc("Automobiles");
	    queue = new LinkedList();
	}
	
	/** Empty constructor for loading */
	public Drivelane() { }
	



















	
	/*============================================*/
	/* Basic GET and SET methods                  */
	/*============================================*/


	/** Returns the Id of this drivelane. */
	public int getId() { return Id; }
	/** Sets the Id of this drivelane. */
	public void setId(int newid) { Id = newid; }
	
	/** Returns the queue of this Drivelane */
	public LinkedList getQueue() { return queue; }
	/** Sets the queue of this Drivelane */
	public void setQueue(LinkedList q) { queue = q; }
	
	/** Returns the type of this Drivelane */
	public int getType() { return type; }
	/** Sets the type of this Drivelane */
	public void setType(int t) { type = t; }

	/** Returns the Road that this Drivelane belongs to */
	public Road getRoad() { return road; }
	/** Sets the Road that this Drivelane belongs to */
	public void setRoad(Road r) { road = r; }

	/** Returns the Sign that regulates the traffic on this Drivelane */
	public Sign getSign() { return sign; }
	/** Sets the Sign that regulates the traffic on this Drivelane */
	public void setSign(Sign s) { sign = s; }

	/** Returns the roads users can move to when crossing the Node this lane leads to */
	public boolean[] getTargets() { return targets; }
	/** Sets the roads users can move to when crossing the Node this lane leads to */
	public void setTargets(boolean[] r) { targets = r; }
	
	/** Return the last cycle this Drivelane has moved its Roadusers */
	public int getCycleMoved() { return cycleMoved; }
	/** Sets the last cycle this Drivelane has moved its Roadusers */
	public void setCycleMoved(int cycle) { cycleMoved = cycle; }

	/** Return the last cycle this Drivelane was last asked about the movements of its Roadusers */
	public int getCycleAsked() { return cycleAsked; }
	/** Sets the last cycle this Drivelane was last asked about the movements of its Roadusers */
	public void setCycleAsked(int cycle) { cycleAsked = cycle; }

	/** Returns the length of the tail of this Drivelane */
	public int getTailLength() {
		if(getNodeComesFrom() instanceof Junction)
			return getNodeComesFrom().getWidth() ; // Nice and easy
		else
			return 0;
	}
	
	/** Returns the name of this drivelane. It is unique, and only used by the GUI for decoration. */
	public String getName() { return "Drivelane " + Id; }

	/** Returns the length of this Drivelane including tail*/
	public int getCompleteLength() {
		return getLength()+getTailLength();
	}



    public int getFreeUnitsInFront(Roaduser ru) {
        int posroom = ru.getPosition();
		ListIterator li = queue.listIterator();
		Roaduser r = null;
		while (li.hasNext()) {
			r = (Roaduser)li.next();
			if(r==ru)
			    break;
			posroom -= r.getLength();
		}
		return posroom;
    }












	/*============================================*/
	/* Selectable                                 */
	/*============================================*/

	public Rectangle getBounds() { return getComplexBounds().getBounds(); }
	public Shape getComplexBounds() {
		Area a = new Area();
		if (bounds != null)
			for (int i=0; i < bounds.length; i++)
				a.add(new Area(bounds[i]));
		return a;
	}
	public int getDistance(Point p) { return (int)getCenterPoint().distance(p); }
	public Point getSelectionPoint() { return getCenterPoint(); }
	public Point[] getSelectionPoints() { return null; }
	public Point getCenterPoint() {
		Rectangle r = getBounds();
		return new Point(r.x + r.width / 2, r.y + r.height / 2);
	}
	public boolean isSelectable() { return true; }
	public boolean hasChildren() { return false; }
	public Enumeration getChildren() { return new ListEnumeration(queue); }
















	/*============================================*/
	/* LOAD and SAVE                              */
	/*============================================*/

	public void load (XMLElement myElement,XMLLoader loader) throws XMLTreeException,IOException,XMLInvalidInputException
	{
		Id = myElement.getAttribute("id").getIntValue();
		cycleMoved = myElement.getAttribute("cycle-moved").getIntValue();
		cycleAsked = myElement.getAttribute("cycle-asked").getIntValue();
		type=myElement.getAttribute("type").getIntValue();
		loadData.roadId = myElement.getAttribute("road-id").getIntValue();
		targets = (boolean[])XMLArray.loadArray(this,loader);
		
		if (loader.getNextElementName().equals("sign-tl"))
			sign=new TrafficLight();
		else if (loader.getNextElementName().equals("sign-no"))
			sign=new NoSign();
		else
			throw new XMLInvalidInputException
     			("A drivelane in road "+loadData.roadId+
			" couldn't load its sign. No sign element found.");
   		loader.load(this,sign);
   		sign.setLane(this);
   		queue=(LinkedList)XMLArray.loadArray(this,loader);
   	}

	public XMLElement saveSelf () throws XMLCannotSaveException
	{
		XMLElement result = new XMLElement("lane");
		result.addAttribute(new XMLAttribute("id", Id));
		result.addAttribute(new XMLAttribute("cycle-moved", cycleMoved));
		result.addAttribute(new XMLAttribute("cycle-asked", cycleAsked));
		result.addAttribute(new XMLAttribute("road-id", road.getId()));
		result.addAttribute(new XMLAttribute("type", type));
		return result;
	}
  
 	public void saveChilds (XMLSaver saver) throws XMLTreeException,IOException,XMLCannotSaveException
	{
		XMLArray.saveArray(targets,this,saver,"targets");
		saver.saveObject(sign);
		XMLUtils.setParentName (new ListEnumeration(queue),getXMLName());
		XMLArray.saveArray(queue,this,saver,"queue");
	}
 
 	public String getXMLName ()
	{ 	return parentName+".lane";
	}
	
	public void setParentName (String parentName)
	{	this.parentName=parentName; 
	}
	 
 	class TwoStageLoaderData 
 	{ 	int roadId;
 	}
 
 	public void loadSecondStage (Dictionary dictionaries) throws XMLInvalidInputException,XMLTreeException
 	{	road=(Road)(((Dictionary)(dictionaries.get("road"))).get(new Integer(loadData.roadId)));
		sign.loadSecondStage(dictionaries);
		XMLUtils.loadSecondStage(new ListEnumeration(queue),dictionaries);
	}





























	/*============================================*/
	/* MODIFYING DATA                             */
	/*============================================*/


	/**
	 * Resets this Drivelane.
	 * This will remove all Roadusers on this lane,
	 * reset the cycleMoved and cycleAsked counters,
	 * and reset the sign.
	 * @see Sign#reset()
	 */
	public void reset() {
		//System.out.println("Resetting lane " + Id);
		queue = new LinkedList();
		cycleMoved = -1;
		cycleAsked = -1;
		sign.reset();
	}


	/**
	 * Adds a Roaduser at the end of this lane
	 *
	 * @param ru The roaduser to add
	 * @throws InfraException if the roaduser could not be added
	 */
	public void addRoaduserAtEnd(Roaduser ru) throws InfraException {
		int pos = getCompleteLength() - ru.getLength();
		addRoaduserAtEnd(ru,pos);
	}

	/**
	 * Adds a Roaduser at the end of this lane
	 *
	 * @param ru The roaduser to add
	 * @param pos The position where the roadusers should be added
	 * @throws InfraException if the roaduser could not be added
	 */
	public void addRoaduserAtEnd(Roaduser ru, int pos) throws InfraException {
		if (!queue.isEmpty()) {
			Roaduser last = (Roaduser) queue.getLast();
			if (last.getPosition() + last.getLength() <= pos) {		
				ru.setPosition(pos);
				queue.addLast(ru);
				return;
			}
			else
				throw new InfraException("Position taken.");
		}
		ru.setPosition(pos);
		queue.addLast(ru);
	}
	
	
	
	
	/**
	 * Adds a Roaduser at a given position to the lane
	 *
	 * @param ru The roaduser to add
	 * @param pos The position at which to add the roaduser
	 * @throws InfraException if the position is taken by another roaduser
	 */
	public void addRoaduser(Roaduser ru, int pos) throws InfraException {
		if (!queue.isEmpty()) {		
			ListIterator li = queue.listIterator();
			Roaduser r = null;
			while (li.hasNext()) {
				r = (Roaduser)li.next();
				if (r.getPosition() <= pos && r.getLength() + r.getPosition() > pos)
					throw new InfraException("Position taken");
				if (r.getPosition() >  pos) {
					if (ru.getLength() > r.getPosition() - pos) throw new InfraException("Position taken");
					li.add(ru);
					break;
				}
			}
			if (pos >= r.getPosition() + r.getLength())
				queue.addLast(ru);
		}
		else {
			queue.addLast(ru);
		}
		ru.setPosition(pos);
	}
	
	/**
	 * Removes Roaduser at start of this lane
	 * 
	 * @return The roaduser removed from the queue
	 * @throw InfraException if there are no roadusers on this lane
	 */
	public Roaduser remRoaduserAtStart() throws InfraException {
		if (queue.isEmpty()) throw new InfraException("No roaduser to remove");
		Roaduser ru = (Roaduser)queue.removeFirst();
		//if (ru.getPosition() != 0) throw new InfraException("First Roaduser not at start of lane");
		ru.setPosition(-1);
		return ru;
	}



























	/*============================================*/
	/* COMPLEX GET                                */
	/*============================================*/


	/** Returns an array of primitive Roaduser types that may roam this Drivelane */
	public int[] getTypes() { return Typeutils.getTypes(type); }
	
	/** Returns if a Roaduser of type ruType may use this Drivelane */
	public boolean mayUse(int ruType) { return (type & ruType) == ruType; }

	/** Returns the length of the Road of this Drivelane */
	public int getLength() 
	{ 
		//System.out.println("Drivelane.getlength returns: " + road.getLength());
		return road.getLength(); 
	}
	/** Returns the first Roaduser on this Drivelane */
	public Roaduser getFirstRoaduser() {
		return (Roaduser)queue.getFirst();
	}

	/** Returns the number of Roadusers that are waiting for the Sign of this Drivelane */
	public int getNumRoadusersWaiting() {
		/* old
		ListIterator li = queue.listIterator();
		Roaduser ru = null;
		int pos = 0;
		int ru_pos;
		int count = 0;
		while (li.hasNext()) {
			ru = (Roaduser) li.next();
			ru_pos = ru.getPosition();
			if (ru_pos > pos) return count;
			else if (ru_pos == pos) {
				pos += ru.getLength();
				count++;
			}
		}*/

		ListIterator li = queue.listIterator();
		Roaduser ru = null;
		int pos = 0;
		int ru_pos;
		int count = 0, cnt_step = 0;
		while (li.hasNext()) {
			ru = (Roaduser) li.next();
			ru_pos = ru.getPosition();
			
			// was:
			// if(ru_pos > pos) return count;
			// nu: waar ru terecht kan komen, moet nog rekening worden gehouden met inloop vakjes
			if (ru_pos - ru.getSpeed() > pos - cnt_step) return count; // Wont be able to wait.
			else if (ru_pos - ru.getSpeed() <= pos - cnt_step) {
				cnt_step += ru_pos - pos; // The free blocks ahead of ru, if everyone moves on.
				pos = ru_pos + ru.getLength();
				count++;
			}
		}

		return count;
	}
   
	/** Returns if this ru is waiting, if so it will update the waiting position of this RoadUser */
	public boolean updateWaitingPosition(Roaduser current_ru) {
		ListIterator li = queue.listIterator();
		Roaduser ru = null;
		int pos = 0;
		int ru_pos;
		int count = 0, cnt_step = 0;
		while (li.hasNext()) {
			ru = (Roaduser) li.next();
      if(ru.equals(current_ru)) {
		  	ru_pos = ru.getPosition();

  			if (ru_pos - ru.getSpeed() > pos - cnt_step) return false; // Wont be able to wait.
  			else if (ru_pos - ru.getSpeed() <= pos - cnt_step) {
  				cnt_step += ru_pos - pos; // The free blocks ahead of ru, if everyone moves on.
  				pos = ru_pos + ru.getLength();
          //Deze wacht dus dus geef zijn huidige wacht positie aan
          /*ru.setLastWaitPointPos(ru_pos);
          ru.setLastWaitPointTl(sign);*/
          return true;
  			}
      }
		}

		return false;
	}


	/** Returns the number of Passengers in the Roadusers that are waiting for the Sign of this Drivelane */
	public int getNumPassengersWaiting() {
		ListIterator li = queue.listIterator();
		Roaduser ru = null;
		int pos = 0;
		int ru_pos;
		int count = 0, cnt_step = 0;
		while (li.hasNext()) {
			ru = (Roaduser) li.next();
			ru_pos = ru.getPosition();
			if (ru_pos - ru.getSpeed() > pos - cnt_step) return count; // Wont be able to wait.
			else if (ru_pos - ru.getSpeed() <= pos - cnt_step) {
				cnt_step += ru_pos - pos; // The free blocks ahead of ru, if everyone moves on.
				pos = ru_pos + ru.getLength();
				count += ru.getNumPassengers();
			}
		}
		return count;
	}
	
	/** Returns the number of blocks taken by Roadusers that are waiting for the Sign of this Drivelane */
	public int getNumBlocksWaiting() {
		ListIterator li = queue.listIterator();
		Roaduser ru = null;
		int pos = 0;
		int ru_pos;
		int count = 0, cnt_step = 0;
		while (li.hasNext()) {
			ru = (Roaduser) li.next();
			ru_pos = ru.getPosition();
			if (ru_pos - ru.getSpeed() > pos - cnt_step) return count; // Wont be able to wait.
			else if (ru_pos - ru.getSpeed() <= pos - cnt_step) {
				cnt_step += ru_pos - pos; // The free blocks ahead of ru, if everyone moves on.
				pos = ru_pos + ru.getLength();
				count += ru.getLength();
			}
		}
		return count;
	}				

	/** Returns the number of blocks taken by Roadusers on this Drivelane */
	public int getNumBlocksTaken() {
		ListIterator li = queue.listIterator();
		Roaduser ru = null;
		int pos = 0;
		int ru_pos;
		int count=0;
		while(li.hasNext()) {
			ru = (Roaduser) li.next();
			count += ru.getLength();
		}
		return count;
	}				

	
	/** Returns the Node that this Drivelane comes from */
	public Node getNodeComesFrom() {
		if (road.getAlphaNode() == getNodeLeadsTo()) return road.getBetaNode();
		else return road.getAlphaNode();
	}

 	/** Returns Node that this Drivelane leads to */
	public Node getNodeLeadsTo() { return sign.getNode(); }

	/** Returns the state of given target */
	public boolean getTarget(int target) throws InfraException {
		if (target < 0 || target > 2) throw new InfraException("Target out of range");
		return targets[target];
	}
	/** Sets the state of given target */
	public void setTarget(int target, boolean state) throws InfraException {
		if (target < 0 || target > 2) throw new InfraException("Target out of range");
		targets[target] = state;
	}

	
	/**
	 * Checks whether length blocks from the given position are free.
	 *
	 * @param position The position in the Queue of this Drivelane.
	 * @param length The amount of blocks that need to be free.
	 */
	public boolean isPosFree(int position, int length) {
		ListIterator li = queue.listIterator();
		Roaduser ru = null;
		int rupos;
		int rulen;
		
		while (li.hasNext()) {
			ru = (Roaduser) li.next();
			rupos = ru.getPosition();
			rulen = ru.getLength();
			
			if (rupos > position + length)
				return true;
			else if (rupos + rulen >= position)
				return false;
		}
		return true;
	}
	
	/**
	 * Checks whether length blocks from the given position, before the given roaduser, are free.
	 *
	 * @param position The position in the Queue of this Drivelane
	 * @param length The amount of blocks that need to be free.
	 * @param me The roaduser
	 * @return whether or not the requested position and <code>length</code> blocks are free for the supplied RU
	 */
	public boolean isPosFree(int position, int length, Roaduser me) {
		ListIterator li = queue.listIterator();
		Roaduser ru = null;
		int rupos;
		int rulen;
		
		while (li.hasNext()) {
			ru = (Roaduser) li.next();
			if (ru != me) {
				rupos = ru.getPosition();
				rulen = ru.getLength();
				
				if (rupos + rulen > position) 
					return false;
				else if (rupos > position + length)
					return true;
			}
			else 
				return true;
		}
		return true;		
	}
	
	
	/** Returns whether or not the 'tail' of this Drivelane has any Roadusers on it. */
	public boolean isTailFree() {
		if(!queue.isEmpty()) {
			Roaduser ru = (Roaduser) queue.getLast();
			if(ru.getPosition()<getLength())
				return true;
			else
				return false;
		}
		return true;
	}
	
	public boolean isLastPosFree(int length) {
		int qSize = queue.size();
		int dLength = getLength();
		
		if(qSize > 0) {
			ListIterator li = queue.listIterator();
			Roaduser ru = null;
			int desired_position = getLength()-length;
			int real_pos = getCompleteLength()-length;
			int dlSize = 0;

			while(li.hasNext()) {
				ru = (Roaduser) li.next();
				dlSize += ru.getLength();
				if(ru.getPosition()+ru.getLength()>real_pos)
					return false;
			}
			if(dlSize+length>getLength()) {
				// See if the current inhabitants already fill up the Drivelane
				return false;
			}
			else {
				// The current inhabitants and length fit on the drivelane
				// Furthermore, there is space on the tail (see return false in while{})
				return true;
			}
		}
		return true;
	}
	
	/**
	 * Returns true when this drivelane is completely full
	 * @ author Chaim Z
	 */
	public boolean isFull()
	{
		ListIterator li = queue.listIterator();
		Roaduser ru = null;
		int position=0;
		int rupos;
		int rulen;
		
		// if empty return false
		if (!li.hasNext()) return false;
		
		while (li.hasNext()) 
		{
			ru = (Roaduser) li.next();
			rupos = ru.getPosition();
			
			if (rupos > position) return false;
	
			rulen = ru.getLength();
			position += rulen;
	  	}
	 	return true;
	}
	
	
	/**
	 * Returns the best reacheable position for the supplied Roaduser on the Queue given in the ListIterator
	 * 
	 * @param li The Queue of this Drivelane represented in a ListIterator. li.previous() is the current RoadUser
	 * @param position The position on the Drivelane of the Roaduser
	 * @param length The amount of blocks that have to be free
	 * @param speed_left the number of 'moves' this Roaduser has left in this turn
	 * @param ru The Roaduser to be checked
	 * @return the number of blocks the Roaduser can move ahead
	 */
	public int getPosFree(ListIterator li, int position, int length, int speed_left, Roaduser ru) {
		int best_pos = position;
		int max_pos = position;
		int target_pos = (position - speed_left > 0) ? position-speed_left : 0;

		// Previous should be 'ru'
		Roaduser prv = (Roaduser) li.previous();
			
		if(prv==ru && li.hasPrevious()) {
			prv = (Roaduser) li.previous();
			max_pos = prv.getPosition()+prv.getLength();
			if(max_pos < target_pos)
				best_pos = target_pos;
			else
				best_pos = max_pos;
			li.next();
		}
		else
			best_pos = target_pos;
		
		li.next();	// Setting the ListIterator back in the position we got it like.

		if(best_pos != position)		// The Roaduser can advance some positions
			return best_pos;
		else
			return 0;
	}
	
	public LinkedList getCompleteQueue()
	{	return queue;
	}							


	
	/** Clears the bounds of this drivelane */
	public void clearCurveBounds()
	{
		bounds = null;
	}

	/** Adds new bounds to the bounds of this drivelane */
	public void addCurveBounds(Shape s)
	{
		if (bounds == null) {
			bounds = new Shape[1];
			bounds[0] = s;
		}
		else
			bounds = (Shape[])Arrayutils.addElement(bounds, s);
	}

	/* Returns true if the given point is contained in this drivelane */
	public boolean contains(Point p)
	{
		if (bounds != null)
			for (int i = 0; i < bounds.length; i++)
			{
				if (bounds[i].contains(p))
					return true;
			}
		return false;
	}
}