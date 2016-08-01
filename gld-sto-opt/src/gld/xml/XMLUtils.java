
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

package gld.xml;

import gld.utils.StringUtils;
import java.util.Dictionary;
import java.util.Enumeration;
import java.util.StringTokenizer;

/** This class contains a few static utility methods that are
  * used by the various parts of the XML parser
 */

public class XMLUtils
{ 
  /** Calculates the last name of a XML path name
    * @param fullName The path name
    * @return Its last name
   */
  public static String getLastName (String fullName)
  { return fullName.substring(fullName.lastIndexOf('.')+1,fullName.length());
  }
  
  /** Calculates the last name of a XMLSerializable object
    * @param object The XMLSerializable object
    * @return Its last name
   */
  public static String getLastName (XMLSerializable object)
  { return getLastName(object.getXMLName());
  }
  
  
  /** Calculates the parent path name of a XML path name
    * @param fullName The path name
    * @return The path of its parent
   */
  public static String getParentName (String fullName)
  { int lastDotPos=fullName.lastIndexOf('.');
    if (lastDotPos==-1)
    { return "";
    }
    else
    { return fullName.substring (0,lastDotPos);
    }
  }
  
  /** Calculates the parent path name of a XML path name
    * @param object The XMLSerializable object
    * @return The path of its parent
   */
  public static String getParentName (XMLSerializable object)
  { return getParentName(object.getXMLName());
  }
  
  
  /** Calculates the level of a XML path name
    * @param fullName The path name
    * @return the level
   */
  public static int getLevel (String fullName)
  { return (new StringTokenizer(fullName,".")).countTokens();
  }

/** Encodes characters which can confuse the XML parser. Every 
    * non-alphanumeric string should be processed by this method before saving!
    * @param s The string to encode
    * @return The result
   */
  public static String toXMLString (String s)
  { if (s==null)
       return ("&null;");
    String string=new String(s);
    char[] from={'&',' ','"','<','>','/'};
    String[] to={"&amp;","&#32;","&quot;","&lt;","&gt;","&#47;"};
    return StringUtils.replaceList(string,from,to);
  }
  
  /** Decodes a string which is encoded with toXMLString. 
    * @param s The string to decode
    * @return The result
   */
  public static String fromXMLString (String s)
  { if (s.equals("&null;"))
       return null;
    String string=new String(s);
    String[] from={"&#47;","&lt;","&gt;","&quot;","&#32;","&#amp"},
             to={"/","<",">","\""," ","&"};
    return StringUtils.replaceList(string,from,to);
  } 
  
  /** Convert a full XML name to its generic name. This is done by checking
    * if each part of the name is complex (if it consists of a generic name plus
    * a specific name). If so, then the specific part is removed. Array elements
    * are also removed from the name.
    * Example of a complex name : model.infrastructure.node-edge.spdata
    * Result of getGenericName  : model.infrastructure.node.spdata
    * @param name The XML name
    * @return The result
   */
  public static String getGenericName (String fullName)
  { 	return getGenericName (fullName,true);
  }
  
  /** Convert a full XML name to its generic name. This is done by checking
    * if each part of the name is complex (if it consists of a generic name plus
    * a specific name). If so, then the specific part is removed. 
    * Example of a complex name : model.infrastructure.node-edge.spdata
    * Result of getGenericName  : model.infrastructure.node.spdata
    * @param name The XML name
	 * @param removeArrays Wether array elements have to be removed from
	 *        the name.
    * @return The result
   */
  public static String getGenericName (String fullName,boolean removeArrays)
  { StringTokenizer s=new StringTokenizer (fullName,".");
    String result="",tmp;
    int dashIndex;
    while (s.hasMoreTokens())
    { tmp=s.nextToken();
      if ( removeArrays && (("array").equals(tmp) || ("element").equals(tmp)))
         continue;
      if (result.length()>0)
         result+=".";
      if ( (dashIndex=tmp.indexOf('-')) == -1)
      { result+=tmp;
      }
      else
      { result+=tmp.substring(0,dashIndex);
      }
    }
    return result;
  }
  
  
  /** Takes the last name of a full XML tag name. If it consists of a generic
    * part and a specific part, then the generic part is removed.
    * examples : "model.infrastructure.node-edge" becomes "edge"
    *         :  "model.infrastructure' becomes "infrastructure"
    * @param fullName The full XML tag name to convert
    * @return The result
   */
  public static String getSpecificLastName (String fullName)
  { int beginIndex=Math.max(Math.max(fullName.lastIndexOf('.'),
                            fullName.lastIndexOf('-')),0)+1;
    return fullName.substring(beginIndex,fullName.length());
  }
  
	/** Load the second stage of an enum_vareration of objects that implement
	  * TwoStageLoader
	  * @param e The enum_vareration of TwoStageLoaders
	  * @param dictionaries The main dictionary to load them with
	  * @throws ClassCastException If one of the elements in the enum_vareration
	  *         is not a TwoStageLoader
	 */
  	public static void loadSecondStage (Enumeration e,Dictionary dictionaries) throws XMLInvalidInputException,XMLTreeException
	{ while (e.hasMoreElements())
		((TwoStageLoader)(e.nextElement())).loadSecondStage(dictionaries);
	}
		
	/** Sets new parentNames for every element in an enum_vareration of
	  * XMLSerializables.
	  * @param e The enum_vareration of XMLSerializables
	  * @param newParentName The new p
	  * @throws ClassCastException If one of the elements in the enum_vareration
	  *         is not XMLSerializable
	 */
  	public static void setParentName(Enumeration e,String newParentName) throws XMLTreeException
	{ while (e.hasMoreElements())
		((XMLSerializable)(e.nextElement())).setParentName(newParentName);
	}
  
}
  