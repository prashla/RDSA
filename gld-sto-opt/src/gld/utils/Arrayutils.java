
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

import java.lang.reflect.Array;
import java.util.Enumeration;
import java.util.NoSuchElementException;
import java.util.Random;

/**
 *
 * These functions can be used to add and remove elements from an array.
 *
 * @author Joep Moritz && Jilles V
 * @version 1.0
 */

public class Arrayutils
{
	/**
	 * Adds an element to the end of an array.
	 * Although the return type is an Object,
	 * you can cast it to an array of the same type as ar.
	 *
	 * @param ar The array to add an element to.
	 * @param elem The element to add to the end of the array.
	 * @return A new array with the same elements as ar, with elem at the end.
	 */
	public static Object addElement(Object ar, Object elem) {
		int arl = Array.getLength(ar);
		Object newar = Array.newInstance(ar.getClass().getComponentType(), arl+1);
		System.arraycopy(ar, 0, newar, 0, arl);
		Array.set(newar, arl, elem);
		return newar;
	}
	
	/**
	 * Adds an element to the end of an array, if and only if
	 * the element is not already in the array. Comparison is done
	 * using the == operator, NOT the Object.equals method.
	 * Although the return type is an Object,
	 * you can cast it to an array of the same type as ar.
	 *
	 * @param ar The array to add an element to.
	 * @param elem The element to add to the end of the array.
	 * @return A new array with the same elements as ar, with elem at the end.
	 */
	public static Object addElementUnique(Object ar, Object elem) {
		if (findElement(ar, elem) == -1)
			return addElement(ar, elem);
		return ar;
	}	
	
	/**
	 * Creates a new array containing the elements of two arrays.
	 * The elements of the first array will be located at the front of the resulting array,
	 * while the elements of the second will appear at the end. This order is guaranteed.
	 * Although the return type is an Object,
	 * you can cast it to an array of the same type as the first array.
	 *
	 * @param ar1 The array to add the other array to.
	 * @param ar2 The array to add at the end of the first array.
	 * @return A new array with the same elements as ar, with all elements of arrtoadd added at the end.
	 */
	public static Object addArray(Object ar1, Object ar2) {
		int arl1 = Array.getLength(ar1);
		int arl2 = Array.getLength(ar2);
		Object newar = Array.newInstance(ar1.getClass().getComponentType(), arl1 + arl2);
		System.arraycopy(ar1, 0, newar, 0, arl1);
		System.arraycopy(ar2, 0, newar, arl1, arl2);
		return newar;
	}

	/**
	 * Creates a new array containing the unique elements of two arrays.
	 * Comparison is done using the == operator, NOT Object.equals.
	 * This means the null element will appear once in the result, iff
	 * it appears at least once in either of the two input arrays.
	 * The order of the elements is the same as the two input arrays, with
	 * the elements of the first array up front.
	 * Although the return type is an Object, you can cast it to an array
	 * of the same type as the first (or second for that matter) array.
	 *
	 * @param ar1 The array to add the other array to.
	 * @param ar2 The array to add at the end of the first array.
	 * @return A new array with the unique elements of the first and second array.
	 */
	public static Object addArrayUnique(Object ar1, Object ar2) {
		int arl1 = Array.getLength(ar1);
		int arl2 = Array.getLength(ar2);
		int pointer = 0;
		Object elem = null;
		Object[] newar = (Object[])Array.newInstance(ar1.getClass().getComponentType(), arl1 + arl2);
		for (int i=0; i < arl1; i++) {
			elem = Array.get(ar1, i);
			if (findElement(newar, elem) == -1) newar[pointer++] = elem;
		}
		for (int i=0; i < arl2; i++) {
			elem = Array.get(ar2, i);
			if (findElement(newar, elem) == -1) newar[pointer++] = elem;
		}
		return cropArray(newar, pointer);
	}
	
	
	/**
	 * Removes the element at given position, and moves all elements above that 1 position down.
	 * Although the return type is an Object[],
	 * you can cast it to an array of the same type as ar.
	 *
	 * @param ar The array to remove the element from.
	 * @param pos The position of the element to remove.
	 * @return A new array containing all elements of ar, except the element at position pos.
	 */
	public static Object remElement(Object ar, int pos) {
		int arl = Array.getLength(ar);
		Object newar = Array.newInstance(ar.getClass().getComponentType(), arl - 1);
		System.arraycopy(ar, 0, newar, 0, pos);
		if (pos != arl - 1) System.arraycopy(ar, pos+1, newar, pos, arl - pos - 1);
		return newar;
	}
	
	/**
	 * Removes an element from an array. Not an equals method is used,
	 * but a direct == comparison. This method will thus only remove the exact same
	 * element from the array.
	 * Although the return type is an Object[],
	 * you can cast it to an array of the same type as ar.
	 *
	 * @param ar The array to remove the element from.
	 * @param elem The element to remove.
	 * @return A new array containing all elements of ar, except elem.
	 */
	public static Object remElement(Object ar, Object elem) {
		int index = findElement(ar, elem);
		if (index == -1) return ar;
		return remElement(ar, index);
	}
	
	/**
	 * Replaces an element in an array with another.
	 * Although the return type is an Object[],
	 * you can cast it to an array of the same type as ar.
	 *
	 * @param ar The array to replace the element in.
	 * @param oldo The element to replace.
	 * @param newo The element to replace it with.
	 * @return A new array containing all elements of ar, except elem.
	 */
	public static Object setElement(Object ar, Object oldo, Object newo) {
		int index = findElement(ar, oldo);
		if (index == -1) return addElement(ar, newo);
		Array.set(ar, index, newo);
		return ar;
	}
		
	
	
	/**
	 * Checks if an element is part of an array, using the == operator.
	 *
	 * @param ar The array to search.
	 * @param elem The element to search for.
	 * @return The index of the element in the array. If the element is not in the array, returns -1.
	 */
	public static int findElement(Object ar, Object elem) {
		int arl = Array.getLength(ar);
		for (int i=0; i < arl; i++)
			if (Array.get(ar, i) == elem) return i;
		return -1;
	}
	
	public static int findElementA(Object[] ar, Object elem) {
		int arl = ar.length;
		for (int i=0; i < arl; i++)
			if (ar[i] == elem) return i;
		return -1;
	}
	
	/*Shuffle all elements in this array */
	public static int[] randomizeIntArray(int[] arr, Random r)
	{
		int temp, switchVal;
		int arlen = arr.length;
		
		for (int i=0; i<arlen; i++) {
			temp = (int) Math.floor(r.nextFloat()*( arlen-1 )) ;
			switchVal = arr[i] ;
			arr[i] = arr[temp] ;
			arr[temp] = switchVal;
		}
		return arr ;
	}

	
	

	/**
	 * Replaces the elements of the provided array in a random fashion.
	 *
	 * @param ar The array to randomize.
	 * @return The array in which the elements are replaced randomly..
	 */
	public static Object[] randomizeArray(Object[] ar) {
		Random r = new Random(1015);
	    
		int newplace, place;
		int numswaps = ar.length;
		Object temp;
		for(;numswaps>0;numswaps--) {
			newplace = (int) Math.round(r.nextFloat()*(ar.length-1));
			place = (int) Math.round(r.nextFloat()*(ar.length-1));
			
			temp = ar[newplace];
			ar[newplace] = ar[place];
			ar[place] = temp;
		}
		return ar;
	}

	/** Change the size of an array. If the size of the array increases,
	  * then fill the new elements of the array (at the end) with new instances of a
	  * specific class.
	  * @param input The original array
	  * @param newLength The new length of the array
	  * @param newClass The class of the new elements. It is necessary that
	  *        objects of this class can be instantiated with
	  *	   Class.newInstance. If the new length is less than the
	  *	   old length of the array, then this parameter is irrelevant.
	  * @returns The result array
	  * @throws ClassNotFoundException If this method wants to create 
	  *         objects, but it cannot find the class that newClass is
	  *         referring to.
	  * @throws InstantiationException If this method wants to create
	  *         objects, but something goes wrong.
	 */
	public static Object adjustArraySize (Object[] input, int newLength,
						Class newClass)
			throws ClassNotFoundException,InstantiationException,IllegalAccessException			
	{ 	if (input.length==newLength)
	  		return input;
		Object result=Array.newInstance(newClass,newLength);
		System.arraycopy(input,0,result,0,
			Math.min(newLength,input.length));
		if ( newLength > input.length)
		{ 	for (int t=input.length;t<Array.getLength(result);t++)
				Array.set(result,t,newClass.newInstance());
		}
		return result;				
	}		

	
	/**
	 * Crops an array to specified length.
	 * Removes all elements after the specified length.
	 *
	 * @param ar The array to crop.
	 * @param len The length to crop to.
	 * @return A new array with the same elements as ar, cropped to length len
	 */
	public static Object cropArray(Object ar, int length) {
		Object newar = Array.newInstance(ar.getClass().getComponentType(), length);
		System.arraycopy(ar, 0, newar, 0, length);
		return newar;
	}
	
	/**
	 * Resizes an array to specified length.
	 *
	 * @param ar The array to resize.
	 * @param len The length to resize to.
	 * @return A new array with the same elements as ar, resized to length len
	 */
	public static Object resizeArray(Object[] ar, int length) {
		Object newar = Array.newInstance(ar.getClass().getComponentType(), length);
		int new_length = length>ar.length?length:ar.length;
		System.arraycopy(ar, 0, newar, 0, new_length);
		return newar;
	}



	/**
	 * Concatenates a two dimensional array to a one dimensional
	 *
	 * @param ar Two dimensional array to concatenate
	 * @return A new array with all the elements of the original array in its second dimension
	 */
	public static Object concatArray(Object ar) {
		int len1 = Array.getLength(ar);
		int len2 = 0;
		int cnt2 = 0;
		
		for (int i=0; i < len1; i++)
			len2 += Array.getLength(Array.get(ar, i));
		
		Object newar = Array.newInstance(ar.getClass().getComponentType().getComponentType(), len2);
		Object temp = null;
		int len = 0;
		
		for (int i=0; i < len1; i++) {
			temp = Array.get(ar, i);
			len = Array.getLength(temp);
			System.arraycopy(temp, 0, newar, cnt2, len);
			cnt2 += len;
		}

		return newar;
	}
	
//	/** Makes a new array of the specified type and length */
/*	public static Object getNewArray (String className,int length) throws ClassNotFoundException
	{ if ( getDimensionClassName(className) == 1)
	  { return Array.newInstance (getComponentClass(className),length);
	  }
	  else
	  { throw new IllegalArgumentException 
	    ("Arrayutils.getNewArray cannot produce new one dimensional array:"+
	     "dimension doesn't match method.");
	  }
	}
*/

//	/** Makes a new array of the specified type,length and width*/
/*
	public static Object getNewArray (String className,int length,int width) throws ClassNotFoundException
	{ if ( getDimensionClassName(className) == 2)
	  { int[] dimension={length,width};
	    return Array.newInstance (getComponentClass(className),dimension);
	  }
	  else
	  { throw new IllegalArgumentException 
	    ("Arrayutils.getNewArray cannot produce new two dimensional array:"+
	     "dimension doesn't match method.");
	  }
	}
*/	
	/** Determines the dimension of an array by the name of its class
	  * @param The classname of an array.
	  * @returns Its dimension
	 */
	public static int getDimensionClassName (String className) 
	{ if (className==null)
	  { return 0;
	  }
	  else 
	  { int index=0,counter=0;
	    while (index<className.length() && className.charAt(index++)=='[')
	          counter++;
	    return counter;
	  }
	}
	
  /** Calculates the *direct* component of an array by the name of its class.
    * So for instance : the direct component of "[[[I" is "[[I".
    * @param The classname of an array
    * @returns The classname of its direct component. If the direct component
    *          is a "normal class" then the classname is converted to the
    *          standard notation. Names of primitive types are *NOT* converted
    *          (see java.lang.Class documentation for their meanings).
   */
  public static String getComponentClassName (String className)
  { String result;
    if (className.startsWith("["))
    { result=className.substring(1,className.length());
    }
    else
    { result=className; // Not an array
    }
    if (result.startsWith ("L") && result.endsWith (";"))
    { return result.substring(1,result.length()-1);
    }
    else
    { return result;
    }
  }

  /** Calculates the *atomary* component of an array by the name of its class.
    * So for instance : the atomary component of "[[[I" is "I".
    * @param The classname of an array
    * @returns The classname of its atomary component. If the direct component
    *          is a "normal class" then the classname is converted to the
    *          standard notation. Names of primitive types are *NOT* converted
    *          (see java.lang.Class documentation for their meanings).
   */
  public static String getAtomaryComponent (String className)
  { String result=className.substring(className.lastIndexOf('[')+1,className.length());
    if (result.startsWith ("L") && result.endsWith (";"))
    { return result.substring(1,result.length()-1);
    }
    else
    { return result;
    }
  }
	
  /** Creates a one-dimensional array of a certain class.
    * @param className The classname of the components of this array. This
    *                  method accepts normal class names and the one letter
    *                  designations for primitive types as used in 
    *                  java.lang.Class.
    * @param length The length of the array. Positive int of course.
    * @returns The resulting array
    * @throws ClassNotFoundException It's all in the name
   */
  public static Object createSingleArray (String className, int length) throws ClassNotFoundException
  { // Yes, this is ugly. 
    if ("B".equals(className))
    { return new byte[length];
    }
    else if ("C".equals(className))
    { return new char[length];
    }
    else if ("D".equals(className))
    { return new double[length];
    }
    else if ("F".equals(className))
    { return new float[length];
    }
    else if ("I".equals(className))
    { return new int[length];
    }
    else if ("J".equals(className))
    { return new long[length];
    }
    else if ("S".equals(className))
    { return new short[length];
    }
    else if ("Z".equals(className))
    { return new boolean[length];
    }
    else
    { return Array.newInstance(Class.forName(className),length);
    }
  }

  /** Creates a n-dimensional array of a certain class
    * @param className The classname of the components of this array. This
    *                  method accepts normal class names and the one letter
    *                  designations for primitive types as used in 
    *                  java.lang.Class.
    * @param length The length of the array. Positive int of course.
    * @param dimension The dimension of this array
    * @throws ClassNotFoundException It's all in the name
    * @returns The resulting array.
   */
  public static Object createArray (String className,int length,int dimension) throws ClassNotFoundException
  { if (dimension<1)
    { return null;
    }
    else if (dimension==1)
    { return createSingleArray(className,length);
    }
    else
    { // Woej! Recursion !
      return  Array.newInstance(createArray(className,length,dimension-1).getClass(),length);
    }
  }
  
  /** Yet another way to create an array
    * @param arrayName The classname of the array that you want to create
    * @param length The desired length
    * @throws ClassNotFoundException It's all in the name
    * @returns The resulting array.
   */
  public static Object createArray (String arrayName,int length) throws ClassNotFoundException
  { return createArray(getAtomaryComponent(arrayName),length,
                       getDimensionClassName(arrayName));
  }
           
                       
	
	/** Returns an Enumeration over the given array */
	public static Enumeration getEnumeration(Object[] ar) {
		return new ArrayEnumerator(ar);
	}
	/** Returns an Enumeration over the given 2d array */
	public static Enumeration getEnumeration(Object[][] ar) {
		return new ArrayEnumerator(ar);
	}
	/** Returns an empty Enumeration */
	public static Enumeration getEmptyEnumeration() {
		return new Enumeration() {
			public boolean hasMoreElements() { return false; }
			public Object nextElement() throws NoSuchElementException {
				throw new NoSuchElementException();
			}
		};
	}

	
	private static class ArrayEnumerator implements Enumeration
	{
		Object[][] ar;
		int i, j;
		
		public ArrayEnumerator(Object[] _ar) {
			ar = new Object[1][];
			ar[0] = _ar;
			i = 0;
			j = 0;
		}
		
		public ArrayEnumerator(Object[][] _ar) {
			ar = _ar;
			i = 0;
			j = 0;
		}
		
		public boolean hasMoreElements() {
			return i < ar.length && j < ar[i].length;
		}
		
		public Object nextElement() throws NoSuchElementException {
			if (!hasMoreElements()) throw new NoSuchElementException();
			Object o = ar[i][j++];
			if (j >= ar[i].length) {
				i++;
				j = 0;
			}
			return o;
		}
	}
}