
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

package gld.algo.edit;

import java.io.*;
import gld.infra.*;
import gld.utils.Typeutils;

/**
 * This class will determine for each node for each exitnode what lanes are
 * on a shortest path towards it.
 */

public class ShortestPathCalculator
{
	private final static int int_MAX = Integer.MAX_VALUE; // To be used as 'infinite'
	private Node[] allNodes;
	private Node[] v;
	private int[] d;
	private final static double derivation_factor = 1.1; // Factor of how much longer a path may be to be seen as a 'shortest path'

	/**
	 * The constructor for a shortest path calculator.
	 */
	public ShortestPathCalculator()
	{
	}

	/**
	 * Calculates all the shortest paths for each EdgeNode in the Infrastructure.
	 * @param infra The Infrastructure of which the shortest paths should be calculated.
	 */

	public void calcAllShortestPaths(Infrastructure infra) throws InfraException
	{	Node[] exits = infra.getSpecialNodes();
		allNodes = infra.getAllNodes();
		int num_exits = exits.length;
		//System.out.println("Calculating Shortest Paths to "+num_exits+" ExitNodes");
		for(int i=0;i<num_exits;i++)
			calcShortestPaths(infra, exits[i]);
	}		

	/**
	 * Calculates the shortest paths to the given Node from each and every other Node in the Infrastructure.
	 * @param infra The Infrastructure of which the shortest paths should be calculated.
	 * @param exit The Node to which the shortest paths should be caluclated.
	 */
	private void calcShortestPaths(Infrastructure infra, Node exit) throws InfraException
	{
		Node node = exit;                                           // The Node variable we need for tracking
		int exitId = exit.getId();                                  // The ID we need for entering the found route in the SPData
		int nearest = int_MAX;                                      // Needed for speeding up the 'findNearest..'
		int nearestId = 0;                                          
		int id = exitId;                                            // The id of the node we're working with
		int num_nodes = infra.getAllNodes().length;                 // Needed for looping through all nodes

		d = new int[num_nodes];										// To keep track of the length of the shortest path, necessary for this algorithm.


		// To get all the types of Roadusers that can reach this ExitNode
        Drivelane[] lanes = node.getInboundLanes();
        Drivelane[] edgelanes = exit.getInboundLanes();
        Drivelane[] splanes;                                        // For looping, no need to recreate the Object all over
        Drivelane l;
        Node newNode;
        int newId;
        int types = 0;
        int lanelength;
        int num_lanes = lanes.length;
        int num_edgelanes = edgelanes.length;
        int num_nodes_left;
        int num_splanes;
        
		for(int i=0;i<num_lanes;i++)
			types |= lanes[i].getType();						// Getting a complete typing
		
		int[] type_arr = Typeutils.getTypes(types);					// Getting each primitive type of that
		int   num_types = type_arr.length;
		
		//System.out.println("--- EdgeNode "+exitId+" ---");
		//System.out.println("Number of inbound lanes on EdgeNode "+exit.getId()+" : "+num_lanes+" totaltype:"+types);

		for(int t=0;t<num_types;t++)								// For each of those Roaduser Types
		{
			
   			for(int i=0;i<num_nodes;i++)						
				d[i] = int_MAX;
			d[exitId] = 0;
			
			lanes = exit.getInboundLanes();
			for(int i=0;i<lanes.length;i++) {
				if(lanes[i].mayUse(type_arr[t])) {
					newNode = lanes[i].getNodeComesFrom();
					newNode.addShortestPath(lanes[i], exitId, type_arr[t],0);
					int tnewId = newNode.getId();
					d[tnewId] = lanes[i].getLength();
					//System.out.println("Adding Shortest Path from "+newNode.getId()+" to:"+exitId+" laneId:"+i+" for type "+type_arr[t]+" : length == "+d[tnewId]);
				}
			}

			
			
			v = (Node[]) allNodes.clone();							// We need all the nodes as we will have to determine the shortest paths
																    // for each node to the given exit.
																    // We use an array as with using the Id as the key, we get an O(1) speed of resolving the actual Node-pointer
			
			num_nodes_left = num_nodes-1;
			v[exitId] = null;
			
			//System.out.println("Now Finding all SPs for Type: "+type_arr[t]);
			
	        while(num_nodes_left != 0)
    	    {
        	    num_nodes_left--;
	            id = findNearestAvailableNode();								// Find the ID of the nearest Node not yet in our 'Cloud'. To help the search, supply a suspect.
	            if (id==-1) break;
    	        node  = v[id];                                          		// Get the pointer to that Node
           	    v[id] = null;	                                         		// Remove that node from our 'todo' list
            	
            	splanes = node.getShortestPaths(exitId,type_arr[t]);			// Get the lanes that are on a shortest path from this node to the exitNode.
            	num_splanes = splanes.length;
            	//System.out.println("At "+id+" to search for SP to "+exitId+" lanes to search: "+num_splanes);

				for(int i=0;i<num_splanes;i++)
				{
					lanes = node.getLanesLeadingTo(splanes[i], type_arr[t]);	// Get an array of Drivelanes that can transfer Roadusers to the shortest-path lane
					int num_working_lanes = lanes.length;
					
					//System.out.println(num_working_lanes+" lanes leading to SPLane "+i);
					
					for(int j=0;j<num_working_lanes;j++)                		// For each of those lanes 
					{
						l = lanes[j];
	                    lanelength = l.getLength();
    	                newNode = l.getNodeComesFrom();
        	            newId   = newNode.getId();
        	            
            	        if(d[id] + lanelength < d[newId])               // Determine if the road to the Exitnode is quicker via node
                	    {                                               // Ifso, Remove the earlier thought-up shortestpaths, and use this one.
                    	    d[newId] = d[id]+lanelength;                // Also update the newly found length of the shortest path.
        	            	//System.out.println("Setted SP from: "+newId+" to "+exitId+" : "+(d[id]+lanelength)+" oldlength: "+d[newId]);
                        	//newNode.setShortestPath(l, exitId, type_arr[t],d[newId]);
                        	
                        	// Instead of 'set' now just add this lane, and remove all 'too long lanes'
							newNode.addShortestPath(l, exitId, type_arr[t], d[newId]);
							newNode.remPaths(exitId, type_arr[t], (int) Math.floor(d[newId]*derivation_factor));
						}
						else if(d[id] + lanelength == d[newId])         // If the length of this found path is the same as one already found,
						{                                               // add this route to the shortest paths.
							//System.out.println("Added SP from: "+newId+" laneId "+j+" to "+exitId+" : "+(d[id]+lanelength)+" oldlength: "+d[newId]);
							newNode.addShortestPath(l, exitId, type_arr[t],d[newId]);                   
                	    }
                	    else if(d[id] + lanelength<=(int)Math.floor(d[newId]*derivation_factor)) {
                	    	newNode.addShortestPath(l, exitId, type_arr[t],d[id]+lanelength);
                	    }
	                }
    	        }
			}
		}
	}
	
	/**
	 *Find the nearest which can be choosen as an alternative destination
	 * @param The Id of the current node
	 * @return The node found
	 */
	public int findNearestAvailableNode(int nodeId)
	{
	    if(v[nodeId]!=null)
	    {
	        return nodeId;
	    }
	    else
	    {
	        //int checked = 0;                      // Possible speedup: Keep track of number of nodes checked in this loop
	        int nearest = int_MAX;                  // continue if less than num_nodes_left;
	        int nearestId = -1;
	        int numNodes = d.length;
	        for(int i=0;i<numNodes;i++)
	        {
	            if(v[i]!=null && d[i] <= nearest)
	            {
	                nearest = d[i];
	                nearestId = i;
	            }
	        }
	        return nearestId;
	    }
	}
	/**
	 *Find the nearest which can be choosen as an alternative destination
	 * @return The node found
	 */
	public int findNearestAvailableNode()
	{
        int nearest = int_MAX;                  // continue if less than num_nodes_left;
        int nearestId = -1;
        int numNodes = d.length;
        for(int i=0;i<numNodes;i++)	{
            if(v[i]!=null && d[i] <= nearest) {
                nearest = d[i];
                nearestId = i;
            }
        }
        ///System.out.println("Nearest is: "+nearestId+" with length:"+nearest);
        return nearestId;
	}


}
