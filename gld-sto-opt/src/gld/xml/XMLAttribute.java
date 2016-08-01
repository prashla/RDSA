
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

/** Container class for attributes of a XML element. Objects have built-in
 *  encapsulation for control-characters.
 */
 
// TODO : - write simple Junit test for encapsulation

public class XMLAttribute
{ protected String name,value;

  /** Make a new XMlAttribute. Encapsulate control characters from the
    * parameters so they don't confuse our XML parser.
    * @param name The name of the attribute
    * @param value The value of the attribute
   */
  public XMLAttribute (String name,String value)
  { this(name,value,true);
  }
  
  /** Make a new XMLAttribute. Choose if you want to encapsulate control chars.
    * @param name The name of the attribute
    * @param value The value of the attribute
    * @param encode If this boolean is true, then control characters in
    *        the parameters will be encapsulated so they don't confuse
    *        our XML parser.
   */
  public XMLAttribute (String name,String value,boolean encode)
  { if (encode)
    { this.name=XMLUtils.toXMLString(name);
      this.value=XMLUtils.toXMLString(value);
    }
    else
    { this.name=name;
      this.value=value;
    }
  }
  
  /** Make a new XMLAttribute. The value is an int 
    * @param name The name of the attribute
    * @param value The value of the attribute
   */
  public XMLAttribute(String name,int value)
  { this (name,value+"",false);
  }

  /** Make a new XMLAttribute. The value is a float
    * @param name The name of the attribute
    * @param value The value of the attribute
   */
  public XMLAttribute(String name,float value)
  { this (name,value+"",false);
  }

  /** Make a new XMLAttribute. The value is a double
    * @param name The name of the attribute
    * @param value The value of the attribute
   */
  public XMLAttribute(String name,double value)
  { this (name,value+"",false);
  }
  
  /** Make a new XMLAttribute. The value is a boolean
    * @param name The name of the attribute
    * @param value The value of the attribute
   */
  public XMLAttribute(String name,boolean value)
  { this(name,(new Boolean(value)).toString(),false);
  }
  
  /** Make a new XMLAttribute. The value is a byte
    * @param name The name of the attribute
    * @param value The value of the attribute
   */
  public XMLAttribute(String name,byte value)
  { this(name,(new Byte(value)).toString(),false);
  }
  
  /** @return A string with name="value". With eventual encapsulation.
   */
  public String toString ()
  { return name+'='+'\"'+value+'\"';
  }
  
  /** @return The name of this attribute
   */
  public String getName ()
  { return XMLUtils.fromXMLString(name);
  }
  
  /** @return The value of this attribute
   */
  public String getValue ()
  { return XMLUtils.fromXMLString(value);
  }
  
  /** @return Tries to convert the value of this attribute to an int and
    *         returns it
    * @throws XMLInvalidInputException If the value of this attribute cannot
    *         be converted to an int.
   */
  public int getIntValue () throws XMLInvalidInputException
  { try
    { return Integer.parseInt(value);
    }
    catch (NumberFormatException e)
    { throw new XMLInvalidInputException
      ("Tried to parse an int in attribute "+name+" , but value "+value+
       "could not be parsed.");
    }
  }
  
  public long getLongValue() throws XMLInvalidInputException
  { try
    { return Long.parseLong(value);
    }
    catch (NumberFormatException e)
    { throw new XMLInvalidInputException
      ("Tried to parse a long in attribute "+name+" , but value "+value+
       "could not be parsed.");
    }
  }

  /** @return Tries to convert the value of this attribute to a boolean and
    *         returns it
    * @throws XMLInvalidInputException If the value of this attribute cannot
    *         be converted to a boolean.
   */
  public boolean getBoolValue () throws XMLInvalidInputException
  { if (value.equals("true"))
       return true;
    if (value.equals("false"))
       return false;
    throw new XMLInvalidInputException   
      ("Tried to parse a boolean in attribute "+name+" , but value "+value+
       "could not be parsed.");
  }

  /** @return Tries to convert the value of this attribute to a float and
    *         returns it
    * @throws XMLInvalidInputException If the value of this attribute cannot
    *         be converted to a float.
   */
  public float getFloatValue () throws XMLInvalidInputException
  { try
    { return Float.valueOf(value).floatValue();
    }
    catch (NumberFormatException e)
    { throw new XMLInvalidInputException   
      ("Tried to parse a float in attribute "+name+" , but value "+value+
       "could not be parsed.");
    }
  }
  
  /** @return Tries to convert the value of this attribute to a double and
    *         returns it
    * @throws XMLInvalidInputException If the value of this attribute cannot
    *         be converted to a double.
   */
  public double getDoubleValue () throws XMLInvalidInputException
  { try
    { return Double.valueOf(value).doubleValue();
    }
    catch (NumberFormatException e)
    { throw new XMLInvalidInputException   
      ("Tried to parse a double in attribute "+name+" , but value "+value+
       "could not be parsed.");
    }
  }
  
  /** @return Tries to convert the value of this attribute to a byte and
    *         returns it
    * @throws XMLInvalidInputException If the value of this attribute cannot
    *         be converted to a byte.
   */
  public byte getByteValue () throws XMLInvalidInputException
  { try
    { return (new Byte(value)).byteValue();
    }
    catch (NumberFormatException e)
    { throw new XMLInvalidInputException   
      ("Tried to parse a byte in attribute "+name+" , but value "+value+
       "could not be parsed.");
    }
  }
  
  
}
