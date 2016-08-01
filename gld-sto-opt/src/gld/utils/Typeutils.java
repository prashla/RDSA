
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

package gld.utils;

/**
 *
 * These functions can be used combine and extract primitive types to and from a type-integer.
 *
 * @author Jilles V
 * @version 1.0
 */

public class Typeutils
{
	/**
	 * Extracts all the primitive types (2^n) from a given type-integer
	 *
	 * @param type The type-integer to extract.
	 * @param returns The array of primitive integer types.
	 */
	public static int[] getTypes(int type) {
		//System.out.println("To Type: "+type);
		if (type == 0) return new int[0];
		int log = (int)Math.floor(Math.log(type) / Math.log(2));
		//System.out.println("Log: "+log);
		int checktype = (int)Math.pow(2, log);
		int[] types = new int[log+1];
		int c = 0;
	
		for (; type > 0; checktype /= 2) {
			//System.out.println("Checktype = "+checktype);
			if (checktype > type) continue;
			types[c] = checktype;
			type -= checktype;
			c++;
		}
		return (int[])Arrayutils.cropArray(types, c);
	}
}