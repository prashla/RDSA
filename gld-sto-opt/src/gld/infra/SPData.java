
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

import gld.utils.Arrayutils;
import gld.xml.*;
import java.io.IOException;
import java.util.Dictionary;
import java.util.Enumeration;
import java.util.Vector;

/**
 *
 * Holds data needed to find the shortest path from
 * a node to an exit node.
 *
 * @author Group Datastructures
 * @version 1.0
 */

public class SPData implements XMLSerializable,TwoStageLoader,InstantiationAssistant
{
	/** The paths known */
	protected Vector paths;
	protected String parentName="model.infrastructure.node";

	public SPData() {
		paths = new Vector(2);
	}
	
	public SPData(int size) {
		paths = new Vector(size);
	}


	/*============================================*/
	/* GETS                                       */
	/*============================================*/
	

	/**
	 * Returns an array of Drivelanes that are on 1 of the shortest paths
	 * from the node this SPData belongs to, to the Node with exiNodeId,
	 * for Roadusers with type ruType.
	 *
	 * @param exitNodeId The Id of the exit node that is your destination.
	 * @param ruType The type of Roaduser.
	 * @return an array of Drivelanes.
	 */
	public Drivelane[] getShortestPaths(int exitNodeId, int ruType) {
		//System.out.println("SPData.Getting shortestPath to:"+exitNodeId+" with type:"+ruType+" from "+paths.size());
		Path p = getPath(exitNodeId, ruType);
		/*System.out.println("SPData.Gotten:"+p.getNodeId()+","+p.getRUType());
		System.out.println("SPData.With "+p.getLanes());*/
		
		if (p != null)
			return p.getLanes();
		else
			return new Drivelane[0];
	}
	
	public int[] getShortestPathDestinations(int ruType) {
		Path p;
		int num_paths = paths.size();
		int counter = 0;
		int[] ps = new int[num_paths];
		for (int i=0; i < num_paths; i++) {
			p = (Path)paths.get(i);
			if (p.getRUType() == ruType) {
				ps[counter] = p.getNodeId();
				counter++;
			}
		}
		if(counter<num_paths)
			return (int[]) Arrayutils.cropArray(ps, counter);
		else
			return ps;
	}

	/*============================================*/
	/* SETS                                       */
	/*============================================*/



	/**
	 * Sets the shortest path to given exit node using roaduser type
	 * to given Drivelane
	 *
	 * @param lane The Drivelane to set the path to
	 * @param exitNodeId The Id of the exit node this path leads to
	 * @param ruType The type of Roaduser
	 */
	public void setShortestPath(Drivelane lane, int exitNodeId, int ruType, int length)
	{
		Path p = getPath(exitNodeId, ruType);
		if(p==null)
			paths.addElement(new Path(exitNodeId, ruType, lane, length));
		else {
			p.empty();
			Drivelane[] lanes = {lane};
			Integer[] lengths = {new Integer(length)};
			p.setLanes(lanes, lengths);
		}
	}
	
	/*============================================*/
	/* ADDS                                       */
	/*============================================*/

  
	/**
	 * Adds a Drivelane to the lanes already found for exitNodeId and ruType.
	 *
	 * @param lane The Drivelane to add to the path
	 * @param exitnodeId The Id of the exit node this path leads to
	 * @param ruType The type of Roaduser
	 */
	public void addShortestPath(Drivelane lane, int exitNodeId, int ruType, int length)
	{
		Path p = getPath(exitNodeId, ruType);
		if(p==null)
			paths.addElement(new Path(exitNodeId, ruType, lane, length));
		else
			p.addLane(lane, length);
	}

	/*============================================*/
	/* REMOVES                                    */
	/*============================================*/

	
	/**
	 * Removes all Drivelanes found for exitNodeId and ruType.
	 *
	 * @param exitnodeId The Id of the exit node this path leads to
	 * @param ruType The type of Roaduser
	 */
	public void remAllPaths(int exitNodeId, int ruType) {
		Path p = getPath(exitNodeId, ruType);
		if (p != null)
			p.empty();
	}
	
	/**
	 * Removes all the Drivelanes found for exitNodeId and ruType and length > length
	 *
	 * @param exitId The Id of the exit node this path leads to
	 * @param ruType The type of Roaduser
	 * @param length The maximum length a path may have to remain
	 */
	public void remPaths(int exitId, int ruType, int length) {
		Path p = getPath(exitId, ruType);
		p.remLanes(length);
	}


	/*============================================*/
	/* PRIVATE                                    */
	/*============================================*/


	/** Gets the Path object for given Node Id and Roaduser type */
	private Path getPath(int exitNodeId, int ruType) {
		Path p;
		for (int i=0; i < paths.size(); i++) {
			p = (Path)paths.get(i);
			if (p.getNodeId() == exitNodeId && p.getRUType() == ruType)
				return p;
		}
		return null;
	}

	/*============================================*/
	/* Load/save                                  */
	/*============================================*/

 	public void load (XMLElement myElement,XMLLoader loader) throws XMLTreeException,IOException,XMLInvalidInputException
 	{ 	paths=(Vector)XMLArray.loadArray(this,loader,this);
 	}

 	public XMLElement saveSelf () throws XMLCannotSaveException
 	{ 	XMLElement result=new XMLElement("spdata");
   	result.addAttribute(new XMLAttribute("num-paths",paths.size()));
  		return result;
 	}
 
 	public void saveChilds (XMLSaver saver) throws XMLTreeException,IOException,XMLCannotSaveException
 	{ 	XMLArray.saveArray(paths,this,saver,"paths");
 	}

 	public String getXMLName ()
 	{ 	return parentName+".spdata";
 	}
 
 	public void setParentName (String parentName)
 	{ 	this.parentName=parentName;
 	}
 
 
 	public void loadSecondStage (Dictionary dictionaries) throws XMLInvalidInputException,XMLTreeException
 	{ 	Enumeration e=paths.elements();
   		while (e.hasMoreElements())
         	((Path)(e.nextElement())).loadSecondStage(dictionaries);
 	}
	
	public boolean canCreateInstance (Class request)
	{	return Path.class.equals(request);
	}
	
	public Object createInstance (Class request) throws 
	      ClassNotFoundException,InstantiationException,IllegalAccessException
	{	if (Path.class.equals(request))
		{	return new Path();
		}
		else
		{ throw new ClassNotFoundException
		  ("SPData IntstantiationAssistant cannot make instances of "+
		   request);
		}
	}
 

	/*============================================*/
	/* Internal classes                           */
	/*============================================*/

  
  /**
   * One Path holds all known Drivelanes that are on a shortest path
   * to a given exitNode using a Roaduser with a certain type,
   * starting at the Node SPData belongs to
   */
	protected class Path implements XMLSerializable,TwoStageLoader
	{
		private int exitNodeId;
		private int min_length = Integer.MAX_VALUE;
		private int max_length = Integer.MAX_VALUE;
		private int ruType;
		private Integer[] lengths;
		private Drivelane[] lanes;
		protected String parentName="model.infrastructure.node.spdata";
		private TwoStageLoaderData loadData=new TwoStageLoaderData();
		
		public Path ()
		{// For loading
		}
    
		Path(int size) {
			lanes = new Drivelane[size];
			lengths = new Integer[size];
		}
		Path(int exitId, int ruT) {
			exitNodeId = exitId;
			ruType = ruT;
			lanes = null;
			lengths = new Integer[0];
		}
		Path(int exitId, int ruT, Drivelane l, int length) {
			exitNodeId = exitId;
			ruType = ruT;
			Drivelane[] nlanes = {l};
			lanes = nlanes;
			Integer[] nlengths = {new Integer(length)};
			lengths = nlengths;
		}
    
		/** Returns all lanes */
		public Drivelane[] getLanes() { return lanes; }
		/** Sets all lanes */
		public void setLanes(Drivelane[] l, Integer[] lens) {
			lanes = l;
			lengths = lens;
		}

		/** Returns the Id of the exitNode */
		public int getNodeId() { return exitNodeId; }
		/** Sets the Id of the exitNode */
		public void setNodeId(int id) { exitNodeId = id; }

		/** Returns the Roaduser type */
		public int getRUType() { return ruType; }
		/** Sets the Roaduser type */
		public void setRUType(int t) { ruType = t; }

		/** Add one Drivelane */
		public void addLane(Drivelane l, int length) {
			int oldlen = lanes.length;
			lanes = (Drivelane[])Arrayutils.addElementUnique(lanes, l);
			if(oldlen<lanes.length) {
				// Something added
				lengths = (Integer[])Arrayutils.addElement(lengths, new Integer(length));
			}
		}
		/** Remove a Drivelane */
		public void remLane(Drivelane l) {
			lanes = (Drivelane[])Arrayutils.remElement(lanes, l);
		}
		
		/** Remove all Drivelanes with pathlength > length */
		public void remLanes(int length) {
			for(int i=0;i<lanes.length;i++) {
				if(lengths[i].intValue()>length) {
					// Removing
					lanes = (Drivelane[])Arrayutils.remElement(lanes, i);
					lengths = (Integer[])Arrayutils.remElement(lengths, i);
				}
			}
		}
		/** Remove all Drivelanes */
		public void empty() {
			lanes = null;
		}
		
		//XMLSerializable implementation
		public void load (XMLElement myElement,XMLLoader loader) throws XMLTreeException,IOException,XMLInvalidInputException
		 { 	ruType=myElement.getAttribute("ru-type").getIntValue();
		   	exitNodeId=myElement.getAttribute("exit-node").getIntValue();
				loadData.laneIds=(int[])XMLArray.loadArray(this,loader);
				lengths=(Integer[])XMLArray.loadArray(this,loader);
		 }

		 public XMLElement saveSelf () throws XMLCannotSaveException
		 { 	XMLElement result=new XMLElement("path");
		   	result.addAttribute(new XMLAttribute("ru-type",ruType));
		   	result.addAttribute(new XMLAttribute("exit-node",exitNodeId));
		   	return result;
		 }
 
		 public void saveChilds (XMLSaver saver) throws XMLTreeException,IOException,XMLCannotSaveException
		 { 	XMLArray.saveArray(getLaneIdArray(),this,saver,"lanes");
		 		XMLArray.saveArray(lengths,this,saver,"lengths");
		 }

		 public String getXMLName ()
		 { 	return parentName+".path";
		 }
		 
		 public void setParentName (String parentName)
		 { 	this.parentName=parentName;
		 }
		 
		 public int[] getLaneIdArray ()
		 { 	int[] result=new int[lanes.length];
		  		for (int t=0;t<lanes.length;t++)
		   	{ if (lanes[t]==null)
		     		result[t]=-1;
		          else
		     		result[t]=lanes[t].getSign().getId();
		     	}
		   	return result;
		 }
		 
		public void loadSecondStage (Dictionary dictionaries) throws XMLInvalidInputException,XMLTreeException
 		{ 	Dictionary laneDictionary=(Dictionary)(dictionaries.get("lane"));
		  	lanes=new Drivelane[loadData.laneIds.length];
		  	for (int t=0;t<loadData.laneIds.length;t++)
		      	lanes[t]=(Drivelane)(laneDictionary.get
		                     (new Integer(loadData.laneIds[t])));		
 		}		 
		
		public class TwoStageLoaderData
		{ 	int[] laneIds;
		}

		
	}
}









/***************
** OLD CODE
**
**
**

	/**
	 * Returns an array of all the Drivelanes that are on a shortest path
	 * from the node this SPData belongs to, to the Node with exiNodeId.
	 *
	 * @param exitNodeId The Id of the exit node that is your destination.
	 * @return an array of Drivelanes.
	 */
/*	public Drivelane[] getShortestPaths(int exitNodeId) {
		
		Vector temp_vector = new Vector();
		Path p;
		int lane_counter = 0;
		for (int i=0; i < paths.size(); i++) {
			p = (Path)paths.get(i);
			if (p.getNodeId() == exitNodeId)
			{
				lane_counter += p.getLanes().length;
				temp_vector.addElement(p);
			}
		}
		int temp_length = temp_vector.size();
		int pos_counter = 0;
		
		// duplicates!
		Drivelane[] lanes = new Drivelane[lane_counter];
		Drivelane[] temp_lanes;
		for (int i=0;i<temp_length;i++)
		{
			p = (Path) temp_vector.elementAt(i);
			temp_lanes = p.getLanes();
			System.arraycopy(temp_lanes, 0, lanes, pos_counter, temp_lanes.length);
			pos_counter += temp_lanes.length;
		}
		temp_lanes  = null;
		temp_vector = null;
		return lanes;
	} */
