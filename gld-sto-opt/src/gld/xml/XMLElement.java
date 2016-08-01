
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
/** This class represents XML elements and their attributes.
 */

import gld.xml.*;
import java.util.*;

// TODO:  - Write JUnit test

public class XMLElement
{ protected String name;
  protected Vector attributes;
  
  /** Make a new XMLElement
    * @param name The name of the new element
   */
  public XMLElement (String name)
  { this.name=name;
    attributes=new Vector();
  }

  /** Make a new XMLElement with the specified name and attributes
    * @param name The name of the new element
    * @param attList An array that contains the attributes for this
    *                element (instances of gld.xml.XMLAttribute)
   */
  public XMLElement (String name,XMLAttribute[] attList)
  { this.name=name;
    attributes=new Vector();
    for (int t=0;t<attList.length;t++)
        attributes.add(attList[t]);
  }
  
  /** @return The name of this XML element
   */
  public String getName ()
  { return name;
  }
  
  /** Change the name of this XML element
    * @param name The new name
   */
  public void setName (String name)
  { this.name=name;
  }
  
  /** Add a new attribute to the attribute list
    * @param attribute The new attribute
   */
  public void addAttribute (XMLAttribute attribute)
  { attributes.add(attribute);
  }
  
  /** Remove the (first) attribute with this name
    * @param name The name of the attribute which has to be removed
    * @throws NoSuchElementException If an attribute with that name
    *         cannot be found.
   */
  public void removeAttribute (String name) throws NoSuchElementException
  { attributes.removeElementAt(getAttributeIndex(name));
  }

  /** Get the index of an attribute in the attribute array
    * @param name The name to search for
    * @throws NoSuchElementIndex If an attribute with that name cannot be found.
    * @return The index of the attribute
   */
  protected int getAttributeIndex (String name) throws NoSuchElementException
  { for (int t=0;t<attributes.size();t++)
    { if (name.equals(((XMLAttribute)(attributes.elementAt(t))).name))
         return t;
    }
    throw new NoSuchElementException 
      ("Could not find XML attribute named "+name);        
  }

  /** Get the (first) attribute with a certain name
   * @param name The name to search for
   * @throws NoSuchElementException If an attribute with that name cannot be found.
   * @return The found attribute
   */
  public XMLAttribute getAttribute (String name) throws NoSuchElementException
  { return (XMLAttribute)(attributes.elementAt(getAttributeIndex(name)));
  }
  
  /* @return A string which describes this XMLElemet
  */
  public String toString ()
  { return "XML element named "+name;
  }

  /* @return An object copy of this XML element
  */
  public Object clone ()
  { return new XMLElement(name,getAttributesArray());
  }
  
  /** Reset the attribute list (make it empty)
   */
  public void removeAllAttributes ()
  { attributes=new Vector ();
  }
  
  /** @return The attribute list in array form
   */
  public XMLAttribute[] getAttributesArray ()
  { XMLAttribute[] result=new XMLAttribute [attributes.size()];
    return (XMLAttribute[])(attributes.toArray(result));
  }
  
  /** Change the attribute with the name of the parameter attribute to the
    * value of the parameter attribute. If an attribute with the name of the
    * parameter doesn't exist, then add the parameter attribute to the
    * attribute list.
    * @param attribute The attribute to set.
   */
  public void setAttribute (XMLAttribute attribute)
  { try
    { ((XMLAttribute)(attributes.elementAt(getAttributeIndex(attribute.name)))).value
                     =attribute.value;
    }
    catch (Exception e)
    { attributes.add(attribute);
    }
  }  
  
  /** @return a string which represents the XML tag which opens this
    *         element in the XML file
   */
  public String getOpenTag ()
  { String result="<"+name;
    XMLAttribute[] attlist=getAttributesArray();
    for (int t=0;t<attlist.length;t++)
        result+=" "+attlist[t].toString();
    result+=">";
    return result;
  }  
  
  /** @return Returns a string which represents the XML tag which closes this
    *         element in the XML file.
   */
  public String getCloseTag ()
  { return "</"+name+">";
  }  
  
  /** Parse a new XML element from a string which contains a XML tag
    * @param string The string to parse
    * @throws XMLTreeException If the string doesn't contain a valid XML tag.
   */
  public static XMLElement parse (String string) throws XMLTreeException
  { StringTokenizer t,a;
    String original=new String(string),attName,attValue; 
    if ( string.startsWith ("</"))
    { throw new XMLTreeException 
      ("XMLElement.parse cannot handle close tags: "+original);
    }
    if (! (string.startsWith("<") && string.endsWith (">")))
    { throw new XMLTreeException 
      ("Invalid XML tag (no < or > ) :"+original);
    }
    string=string.substring(1,string.length()-1);
    t=new StringTokenizer(string," ");
    if (! t.hasMoreTokens() )
    { throw new XMLTreeException
      ("XML tag without valid name "+original);
    }
    XMLElement result=new XMLElement(t.nextToken());
    while (t.hasMoreTokens())
    { a=new StringTokenizer(t.nextToken(),"=");
      if (a.countTokens()!=2)
      { throw new XMLTreeException
        ("Wrong attribute in XML tag :"+original);
      }
      attName=a.nextToken();
      attValue=a.nextToken();
      if (! (attValue.startsWith("\"") && attValue.endsWith ("\"")))
      { throw new XMLTreeException
        ("Wrong attribute value in XML tag (no quotes?) :"+attValue+"/"+original);
      }
      result.addAttribute(new XMLAttribute
             (attName,
	      attValue.substring(1,attValue.length()-1),
	      false));
    }
    return result;	      
  }
  
  
}
