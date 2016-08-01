
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

package gld.algo.dp;
/** This class can be used to create instances of Driving Policies
  * for a specific SimModel and TrafficLightController
 */

import gld.algo.tlc.TLController;
import gld.infra.InfraException;
import gld.infra.Infrastructure; 
import gld.sim.SimModel;
import gld.utils.StringUtils;
import java.util.Dictionary;
import java.util.Hashtable;
import java.util.Random;
 
public class DPFactory
{
	protected SimModel model;
	protected TLController tlc;
	public static final int
		SHORTEST_PATH=0,
		SMARTER_SHORTEST_PATH=1,
		AGGRESSIVE=2,
		COLEARNING=3;
		
	protected static final String[] dpDescs = {
		"Normal shortest path",
		"Least busy shortest path",
		"Aggressive",
		"Colearning"		
		};
	
	protected static final String[] xmlNames = {
		ShortestPathDP.shortXMLName,
		SmarterShortestPathDP.shortXMLName,
		AggressiveDP.shortXMLName,
		ColearnPolicy.shortXMLName		
	};
		
	/** Makes a new DPFactory for a specific SimModel and TLC 
	  * @param model The SimModel to create the algorithm for
	  * @param tlc The traffic light controller to co-operate with
	 */
  	public DPFactory(SimModel model,TLController tlc)
	{ 	this.model=model;
		this.tlc=tlc;
	}
  
	/** Looks up the id of a DP algorithm by its description
	  * @param algoDesc The description of the algorithm
	  * @return The id of the algorithm
	  * @throws NoSuchElementException If there is no algorithm with that
	  *        description.
	 */
	public static int getId (String algoDesc)
	{ 	return StringUtils.getIndexObject(dpDescs,algoDesc);
	}
	
	/** Returns an array of driving policy descriptions */
	public static String[] getDescriptions() { return dpDescs; }
	
	/** Returns a new DrivingPolicy of the requested ID */
	public DrivingPolicy genDP(int dp) throws ClassNotFoundException
	{	return getInstance(dp);
	}
	

  	/** Look up the description of a DP algorithm by its id 
	  * @param algoId The id of the algorithm
	  * @returns The description
	  * @throws NoSuchElementException If there is no algorithm with the
	  *	    specified id.
	*/
  	public static String getDescription (int algoId)
  	{ 	return (String)(StringUtils.lookUpNumber(dpDescs,algoId));
  	}
  
  	/** Gets the number of an algorithm from its XML tag name */
  	public static int getNumberByXMLTagName(String tagName) 
  	{ 	return StringUtils.getIndexObject(xmlNames,tagName);
  	}
  
  	/** Gets a new instance of an algorithm by its number. This method
    	  * is meant to be used for loading.
   	*/
  public DrivingPolicy getInstance (int algoId) throws ClassNotFoundException
  { 	switch (algoId)
    	{ 	case SHORTEST_PATH : return new ShortestPathDP(model,tlc);
			case SMARTER_SHORTEST_PATH : return new ShortestPathDP(model,tlc);
			case AGGRESSIVE : return new AggressiveDP (model,tlc);
			case COLEARNING : return new	ColearnPolicy(model,tlc);
    	}
   	throw new ClassNotFoundException
    		("The DPFactory can't make DP's of type "+algoId);
  }
}  
