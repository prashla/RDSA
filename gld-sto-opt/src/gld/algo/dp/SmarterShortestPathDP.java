
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

import gld.algo.dp.*;
import gld.algo.tlc.*;
import gld.infra.*;
import gld.sim.*;
import gld.xml.*;
import java.io.IOException;
import java.util.Random;

/**
 *
 * This extension of {@see gld.DrivingPolicy} selects the next lane 
 * by finding one which on the shortest path to road user's destination.
 *
 * @author Group Algorithms
 * @version 1.0
 */
public class SmarterShortestPathDP extends DrivingPolicy
{
	public static final String shortXMLName="dp-ssp";
	
	/**
	 * The constructor for a shortest driving policy. 
	 * @param m The model which is used
	 */
    
	public SmarterShortestPathDP(SimModel sim, TLController _tlc) {
		super(sim,_tlc);
	}
    

	/**
	 * The lane to which a car continues his trip.
	 * @param r The road user being asked.
	 * @param allOutgoing All the possible outgoing lanes
	 * @param shortest All the lanes which are in a shortest path to the car's destination
	 * @return The chosen lane.
	 */
	public Drivelane getDirectionLane(Roaduser r, Drivelane lane_now, Drivelane[] allOutgoing, Drivelane[] shortest)
	{	//Create a subset from the 2 sets allOutgoing and shortest
		Drivelane current;
		Drivelane best_lane = null;
		int best_waiting = Integer.MAX_VALUE;
		int num_outgoing = allOutgoing.length;
		int num_shortest = shortest.length;

		for(int i=0; i<allOutgoing.length; i++) {
		    current = allOutgoing[i];
			for(int j=0; j<shortest.length; j++) {
			    if(current.getId() == shortest[j].getId()) // Dirty equallity :( 
				    if(current.getNumRoadusersWaiting()<best_waiting) {
				        best_lane = shortest[j];
					    best_waiting = current.getNumRoadusersWaiting();
				    }	
		    }
	    }
	    return best_lane;
	}	
	
	// Trivial XMLSerializable implementation
	
	public void load (XMLElement myElement,XMLLoader loader) throws XMLTreeException,IOException,XMLInvalidInputException
	{ 	System.out.println("DP SSP loaded");
	}
	
	public XMLElement saveSelf () throws XMLCannotSaveException
	{ 	return new XMLElement(shortXMLName);
	}
  
	public String getXMLName ()
	{ 	return "model."+shortXMLName;
	}
	
}
