
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

import gld.GeneralSettings;
import gld.infra.*;
import java.awt.Color;
import java.util.Random;

/** 
 * This class is used to generate new Roadusers
 *
 * @author Group Datastructures
 * @version 0.01
 */

public class RoaduserFactory
{
	public static float PacChance = 0.0f;
	public static boolean UseCustoms = false;
	protected static Random rnd = new Random(1015);      
	
	/** All the individual types of each concrete Roaduser */
	public static final int CAR=1,
	                        BUS=2,
	                        BICYCLE=4;
	/** All concrete Roaduser types. */
	protected static final int[] concreteTypes = {CAR,BUS,BICYCLE};
	/** *All* roadusers. */
	protected static final int[] ALL = concreteTypes;
	/** All the roadusers that fit in the Automobile category */
	protected static final int[] AUTOMOBILE = {CAR,BUS};
	/** All the roadusers that move on foot */
	protected static final int[] PEDESTRIANS = {};
	/** All the other roadusers that move unmotorized */
	protected static final int[] UNMOTORIZED = {BICYCLE};
	/** All possible Roaduser types. */
	public static final int[] types = { getSuperType(ALL),
																				 getSuperType(AUTOMOBILE),
	                                    getSuperType(UNMOTORIZED),
	                                    CAR,
	                                    BUS,
	                                    BICYCLE}; 
	/** Descriptions of all possible Roaduser types, corresponding to the types in <code>types</code>. */
	public static final String[] typeDescs = {"All roadusers","Automobiles", "Bicycles", "Car", "Bus", "Bicycle"};
	/** Descriptions of all concrete Roaduser types, corresponding to the types in <code>concreteTypes</code>. */
	public static final String[] concreteTypeDescs = {"Car", "Bus", "Bicycle"};
	/** Colors for all concrete Roaduser types, corresponding to the types in <code>concreteTypes</code>. */
	public static final Color[] concreteTypeColors = {Color.blue, Color.red, Color.green};

	private RoaduserFactory() { } // cannot instantiate this


	/** Returns all possible types of roadusers */
	public static int[] getTypes() { return types; }
	/** Returns the number of possible types of roadusers. */
	public static int getNumTypes() { return types.length; }
	/** Returns descriptions of all possible types of roadusers */
	public static String[] getTypeDescs() { return typeDescs; }
	
	/** Returns all concrete types of roadusers. */
	public static int[] getConcreteTypes() { return concreteTypes; }
	/** Returns the number of all concrete types of roadusrs. */
	public static int getNumConcreteTypes() { return concreteTypes.length; }
	/** Returns descriptions of all concrete types of roadusers. */
	public static String[] getConcreteTypeDescs() { return concreteTypeDescs; }
	
	
	/** Generate a new RoadUser
	  * @param type The type number of the Roaduser
	  * @param start The node where the Roadusers should start.
	  * @param dest The destination node of the Roaduser
	  * @param pos The position of the Roaduser in its Drivelane
	  * @throws InfraException If the Roaduser cannot be generated
	  * @returns The generated Roaduser
	 */
 	public static Roaduser genRoaduser(int type, Node start, Node dest, int pos) throws InfraException
	{
		if (UseCustoms)
		{
			Roaduser ru = CustomFactory.genRoaduser(type, start, dest, pos);
			if (ru != null) return ru;
		}
		switch (type) {
			case CAR : return (rnd.nextFloat() < PacChance ? new PacCar(start,dest,0) : new Car(start,dest,0));
			case BUS : return new Bus(start,dest,0);
			case BICYCLE : return new Bicycle(start,dest,0);
		}
		
		throw new InfraException ("The RoaduserFactory couldn't make a Roaduser of type "+type+", reason : unknown class.");
	}
	
	/** Generate a new RoadUser
	  * @param type The type number of the Roaduser
	  * @throws InfraException If the Roaduser cannot be generated
	  * @returns The generated Roaduser
	 */
	public static Roaduser genRoaduser(int type) throws InfraException
	{
		return genRoaduser(type, null, null, 0);
	}

	/** 
	* Seeks the type number of a description of a Roaduser
	* @param desc Description of the Roaduser type.
	* @returns The typenum_varber
	*/
	public static int getTypeByDesc(String desc) {
			for(int i=0; i<types.length; i++)
				if(desc.equals(typeDescs[i]))
					return types[i];
			return 0;
	}

	/** 
	* Returns the desc belonging to the type given.
	* @param type The Roaduser type to return description of.
	*/
	public static String getDescByType(int type) {
			for(int i=0; i<types.length; i++)
				if(type == types[i])
					return typeDescs[i];
			return "";
	}

	/** 
	* Returns the desc belonging to the concrete type given.
	* @param type The Roaduser type to return description of.
	*/
	public static String getDescByConcreteType(int type) {
			for(int i=0; i<concreteTypes.length; i++)
				if(type == concreteTypes[i])
					return concreteTypeDescs[i];
			return "";
	}

	/** 
	* Returns the color belonging to the type given.
	* Colors are only defined for concrete types!
	*/
	public static Color getColorByType(int type) {
		for(int i=0; i<concreteTypes.length; i++)
			if(type == concreteTypes[i])
				return concreteTypeColors[i];
		return null;
	}
	
	protected static int getSuperType(int[] types) {
		int type = 0;
		for(int i=0;i<types.length;i++) {
			type |= types[i];
		}
		return type;
	}

	// Statistical support methods.

	/** 
		* Returns the length a statistics array should have: the number of concrete
		* roadusers plus one (for all roadusers). 
		*/
	public static int statArrayLength() { return concreteTypes.length + 1; }
	
	public static int getStatIndexByType(int ruType) {
		return ruType == 0 ? 0 : (int)(Math.log(ruType)/Math.log(2)+1); 
	}
	
	/** Returns the statIndex belonging to a given concrete (!) ruType. */
	public static int ruTypeToStatIndex(int ruType) {
		return getStatIndexByType(ruType);
	}
	
	/** Returns the concrete ruType belonging to the given statIndex. */
	public static int statIndexToRuType(int index) 
	{ 
		return index == 0 ? 0 : (int)Math.pow(2,index-1);
	}	
	
	/** Returns the roaduser description belonging to the given statIndex. */
	public static String getDescByStatIndex(int index)
	{
		return index == 0 ? "All" : getDescByConcreteType(statIndexToRuType(index));
	}
	
	/** Returns the statIndex belonging to a given roaduser description. */
	public static int getStatIndexByDesc(String desc)
	{
		return ruTypeToStatIndex(getTypeByDesc(desc));
	}
}