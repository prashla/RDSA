
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

import gld.algo.dp.*;
import gld.algo.tlc.*;
import gld.infra.*;
import gld.utils.*;
import java.io.IOException;
import java.lang.reflect.Array;
import java.util.*;

/** Class to convert arrays of arbitrary dimensions from/to XML :
  * Accepted componenttypes :
  * - int,float,double,boolean,byte and their respective encapsulator classes,
  *   Strings and objects that accept XMLSerializable
 */
 
public class XMLArray implements XMLSerializable
{ // The data itself
  protected Object data; 
  protected Class dataClass;
  protected int containerType;
  // Name of parent
  protected String parentName;
  // Short description of what this array contains. Mainly for debugging
  // and error checking.
  protected String alias;
  // InstantiationAssistant for creating objects in XMLArrayElements
  protected InstantiationAssistant watson;
  
  
  // Component type
  protected final static int 
		INT=0, 
		DOUBLE=1,
		FLOAT=2,
		BOOLEAN=3,
		INT_CLASS=4,
		DOUBLE_CLASS=5,
		FLOAT_CLASS=6,
		BOOLEAN_CLASS=7,
		STRING=8,
		XMLSERIALIZABLE=9,
		OTHER_ARRAY=10,
		BYTE=11,
		BYTE_CLASS=12;
		
  protected final static int
  		CONTAINER_ARRAY=0,
		CONTAINER_VECTOR=1,
		CONTAINER_LINKEDLIST=2;		
		
  protected String[] lookUpTable =
  { "I","d","F","Z","java.lang.Integer","java.lang.Double",
    "java.lang.Float","java.lang.Boolean","java.lang.String","L","[","B",
	 "java.lang.Byte"};
    
  /** Make a new XMLArray (constructor for loading)
    * @param parentName The XML name of the object that is going to save
    *        this array.
    * @returns The new XMLArray
   */
  public XMLArray () 
  { 	// Empty
  }
  
  /** Empty constructor for loading 
    * @param parentName The XML name of the object that is going to load
    *        this array.
   *  @returns A new XMLArray which can be used for loading
   */
  public XMLArray (String parentName)
  { this.parentName=parentName;
  }
  
  /** Constructor for loading 
    * @param parentName The XML name of the object that is going to load
    *        this array.
    * @param assistant The IntantiationAssistant to use when creating new
    	     objects (this parameter is only meaningful for XMLArrays that
	     directly or indirectly contain XMLSerializable objects.)
   *  @returns A new XMLArray which can be used for loading
   */
  public XMLArray (String parentName,InstantiationAssistant assistant)
  { this.parentName=parentName;
    watson=assistant;
  }
  

  /** Make a new XMLArray (constructor for saving)
    * @param input The array,Vector or LinkedList to save
    * @param parentName The XML name of the object that is going to save
    *        this array.
    * @throws XMLCannotSaveException If there is something wrong with the input
    *        array.
    * @returns The new XMLArray
   */
  public XMLArray (Object input,String parentName) throws XMLCannotSaveException
  { this(input,parentName,"none"); 
  }
  
  /** Make a new XMLArray (constructor for saving)
    * @param input The array,Vector or LinkedList to save
    * @param parentName The XML name of the object that is going to save
    *        this array.
    * @param alias A short description of what this array contains. Mainly
    *        useful for debugging.
    * @throws XMLCannotSaveException If there is something wrong with the input
    *        array.
    * @returns The new XMLArray
   */
  public XMLArray (Object input,String parentName,String alias) throws XMLCannotSaveException
  { if ( input.getClass().isArray() )
    { containerType=CONTAINER_ARRAY;
    }
    else if ( input instanceof Vector )
    { containerType=CONTAINER_VECTOR;
      input=((Vector)(input)).toArray();
    }
    else if ( input instanceof LinkedList )
    { containerType=CONTAINER_LINKEDLIST;
      input=((LinkedList)(input)).toArray();
    }
    else
    {  	throw new XMLCannotSaveException
    	("Cannot construct XMLArray : input is not of valid type. Input class :"+
	 dataClass.getName());
    }
    data=input;
    dataClass=input.getClass();
    this.parentName=parentName;        
    this.alias=alias;
  }
    
  /** Return the result array,Vector or LinkedList as an object. The XMLArray
    *  has to be loaded first.
    * @throws XMLInvalidInputException If the XMLArray was not loaded or if
    *         there is an internal problem.
    * @returns The array which has been parsed
   */
  public Object getResult () throws XMLInvalidInputException
  { if (data==null)
       throw new XMLInvalidInputException
       ("Tried to get result from XMLArray when it was not loaded in "+
        parentName);
    switch (containerType)
    { case CONTAINER_ARRAY      : return data;
      case CONTAINER_VECTOR     : return getResultVector();
      case CONTAINER_LINKEDLIST : return getResultLinkedList();
    }
   throw new XMLInvalidInputException
   ("Internal problem in XMLArray: result has unknown container ("+
    containerType+")");
  }
  
  /** Gets the result array in Vector form. Used by getResult()
   */
  protected Vector getResultVector ()
  { Object[] array=(Object[])(data);
    Vector result=new Vector(array.length);
    for (int t=0;t<array.length;t++)
        result.add(array[t]);
    return result;
  }
  
  /** Gets the result array in LinkedList form. Used by getResult()
   */
  protected LinkedList getResultLinkedList ()
  { Object[] array=(Object[])(data);
    LinkedList result=new LinkedList();
    for (int t=0;t<array.length;t++)
        result.add(array[t]);
    return result;
  }
  
  /** @returns The instantiation assistant of this XMLArray, or null if doesn't
    *          have one
	*/
  public InstantiationAssistant getInstantiationAssistant ()
  {	return watson;
  }	
  
  /** Sets the InstantiationAssistant of this XMLArray
    * @param The new InstantiationAssistant
	*/
  public void setInstantiationAssistant (InstantiationAssistant assistant)
  {	watson=assistant;
  }	
  
  /** Return the result array,Vector or LinkedList as an object. The XMLArray
    * has to be loaded first. First checks if the result class equals the
    * parameter class.
    * @param expectedClass The class with which the result has to match
    * @throws XMLInvalidInputException If the XMLArray was not loaded, if
    *         there is an internal problem or if the result class doesn't 
    *         match the parameter class.
    * @returns The array which has been parsed
   */
  public Object getResult (Class expectedClass) throws XMLInvalidInputException
  { Object result=getResult();
    if (result.getClass().equals(expectedClass))
    { return result;
    }
    else  
    { throw new XMLInvalidInputException
      ("XMLArray could not return result : expected class "+
        expectedClass.getName()+" didn't match result class "+dataClass.getName());
    }
  }
  
  /** This method tries to guess the component type of the data array. If it 
    * cannot determine the component type, then it assumes that all components
    * are XMLSerializable. This assumption may be wrong, but it doesn't hurt
    * since the saver later checks all components individually.
    * @returns The component type of the data (the result is not 100% reliable)
   */
  protected int getComponentType ()
  { String componentName=Arrayutils.getComponentClassName(dataClass.getName());
    for (int t=0;t< lookUpTable.length;t++)
    { if (componentName.startsWith(lookUpTable[t])) 
         return t;
    }
    if (componentName.startsWith("java.util.Vector"))
    	return XMLArray.OTHER_ARRAY;
    else if (componentName.startsWith("java.util.LinkedList"))
    	return XMLArray.OTHER_ARRAY;
    else	
    	return XMLArray.XMLSERIALIZABLE; // Educated guess
  }  
  
  /** Save an arbitrary array/vector/LinkedList
    * @param array An array to save (Vectors and LinkedLists are also accepted)
    * @param parent The object which is trying to save this array
    * @param saver The XMLSaver to which the array has to be saved
    * @param description A short description of the array
   */
  public static void saveArray (Object array,XMLSerializable parent,
  				XMLSaver saver,String description) 
				throws IOException,XMLCannotSaveException,XMLTreeException
  { saver.saveObject(new XMLArray(array,parent.getXMLName(),description));
  }				
  
  /** Save an arbitrary array/vector/LinkedList
    * @param array An array to save (Vectors and LinkedLists are also accepted)
    * @param parent The object which is trying to save this array
    * @param saver The XMLSaver to which the array has to be saved
   */
  public static void saveArray (Object array,XMLSerializable parent,
  				XMLSaver saver) 
				throws IOException,XMLCannotSaveException,XMLTreeException
  { saveArray(array,parent,saver,"none");
  }				
  
 
  /** Load an arbitrary array. Check if the result matches the expected class.
    * @param parent The object that is loading this array
    * @param loader The XMLLoader to load this array from
    * @param expectedClass The expected class of the result
    * @returns The resulting array,vector or LinkedList
   */
  public static Object loadArray (XMLSerializable parent,XMLLoader loader,
                                  Class expectedClass) 
				  throws XMLInvalidInputException,XMLTreeException,IOException
  { XMLArray tmp=new XMLArray(parent.getXMLName());
    loader.load(parent,tmp);
    return tmp.getResult(expectedClass);
  }
  
  /** Load an arbitrary array. Do not check the result.
    * @param parent The object that is loading this array
    * @param loader The XMLLoader to load this array from
    * @returns The resulting array,vector or LinkedList
   */
  public static Object loadArray (XMLSerializable parent,XMLLoader loader) 
         throws XMLInvalidInputException,XMLTreeException,IOException
  { XMLArray tmp=new XMLArray(parent.getXMLName());
    loader.load(parent,tmp);
    return tmp.getResult();
  }
  
  /** Load an arbitrary array. Instantiate child objects via the 
    * parameter InstantiationAssistant
    * @param parent The object that is loading this array
    * @param loader The XMLLoader to load this array from
    * @param assistant The InstantiationAssistant to use
    * @returns The resulting array,vector or LinkedList
   */
  public static Object loadArray (XMLSerializable parent,XMLLoader loader,
   				  InstantiationAssistant assistant) 
         throws XMLInvalidInputException,XMLTreeException,IOException
  { XMLArray tmp=new XMLArray(parent.getXMLName(),assistant);
    loader.load(parent,tmp);
    return tmp.getResult();
  }
  
  

  // XMLSerializable implementation
  
  public void load (XMLElement myElement,XMLLoader loader) throws XMLTreeException,IOException,XMLInvalidInputException
  { // Create array
    parentName=myElement.getAttribute("parent-name").getValue();
    String className=myElement.getAttribute("class").getValue();
    int componentType=-1;
    containerType=myElement.getAttribute("container").getIntValue();
    int length=myElement.getAttribute("length").getIntValue(); 
    componentType=myElement.getAttribute("component-type").getIntValue();
    try
    {  data=Arrayutils.createArray(className,length);
       dataClass=data.getClass();
    }
    catch (ClassNotFoundException e)
    { throw new XMLInvalidInputException
      ("A new XMLArray couldn't create an instance of its data class:"+
       className+" ,reason :"+e);
    }
    // Populate it
    XMLArrayElement loadObject;
    for (int x=0;x<length;x++)
    { loadObject=new XMLArrayElement();
      loader.load(this,loadObject);
      switch (componentType)
      { case INT: Array.setInt(data,x,loadObject.getIntValue()); break;
        case DOUBLE: Array.setDouble(data,x,loadObject.getDoubleValue()); break;
        case FLOAT: Array.setFloat(data,x,loadObject.getFloatValue()); break;
        case BOOLEAN: Array.setBoolean(data,x,loadObject.getBooleanValue()); break;
        case INT_CLASS: Array.set(data,x,loadObject.getIntClassValue()); break;
        case DOUBLE_CLASS: Array.set(data,x,loadObject.getDoubleClassValue()); break;
        case FLOAT_CLASS: Array.set(data,x,loadObject.getFloatClassValue()); break;
        case BOOLEAN_CLASS: Array.set(data,x,loadObject.getBooleanClassValue());  break;
        case XMLSERIALIZABLE: Array.set(data,x,loadObject.getObjectValue()); break;
        case OTHER_ARRAY: Array.set(data,x,loadObject.getObjectValue()); break;
		  case BYTE: Array.setByte(data,x,loadObject.getByteValue()) ;break;
		  case BYTE_CLASS: Array.set(data,x,loadObject.getByteClassValue());break;
      }
    }		
  }

  public XMLElement saveSelf () throws XMLCannotSaveException
  { XMLElement result=new XMLElement (XMLUtils.getLastName(this));
    result.addAttribute(new XMLAttribute("length",Array.getLength(data)));
    result.addAttribute(new XMLAttribute("class",data.getClass().getName()));
    result.addAttribute(new XMLAttribute("component-type",getComponentType()));
    result.addAttribute(new XMLAttribute("parent-name",parentName));
    result.addAttribute(new XMLAttribute("container",containerType));
    result.addAttribute(new XMLAttribute("name",alias));
    return result;
  }
  
  public void saveChilds (XMLSaver saver) throws XMLTreeException,IOException,XMLCannotSaveException
  { Object toSave;
    int componentType=getComponentType();
    for (int x=0; x< Array.getLength(data); x++)
    {toSave=null;
     switch (componentType)
     { case INT: toSave= (new XMLArrayElement (Array.getInt(data,x))); break;
       case DOUBLE: toSave= (new XMLArrayElement (Array.getDouble(data,x))); break;
       case FLOAT: toSave= (new XMLArrayElement (Array.getFloat(data,x))); break;
       case BOOLEAN: toSave= (new XMLArrayElement (Array.getBoolean(data,x))); break;
       case INT_CLASS: toSave= (new XMLArrayElement ((Integer)Array.get(data,x))); break;
       case DOUBLE_CLASS: toSave= (new XMLArrayElement ((Double)Array.get(data,x))); break;
       case FLOAT_CLASS: toSave= (new XMLArrayElement ((Float)Array.get(data,x))); break;
       case BOOLEAN_CLASS: toSave= (new XMLArrayElement ((Boolean)Array.get(data,x))); break;
       case XMLSERIALIZABLE: toSave= (new XMLArrayElement 
            ((XMLSerializable)((Array.get(data,x))))); break;
       case OTHER_ARRAY: toSave= (new XMLArrayElement
            (new XMLArray (Array.get(data,x),parentName))); break;
		 case BYTE : toSave= (new XMLArrayElement(Array.getByte(data,x)));break;
		 case BYTE_CLASS: toSave= (new XMLArrayElement((Byte)Array.get(data,x)));break;
     }
     if ( toSave==null)
     { System.out.println("WARNING: couldn't save array element in "+
                          getXMLName());
     }
     else
     { saver.saveObject((XMLSerializable)(toSave));
     }
    }
  }
  
  public String getXMLName ()
  { 	return parentName+".array";
  }
  
  public void setParentName (String parentName)
  { 	this.parentName=parentName;
  }
    
  /** Inner class which holds the elements of the array */
  class XMLArrayElement implements XMLSerializable
  { String value="none";
    boolean nullValue=false;
    int type;
    XMLSerializable child=null;
    Object result; // For loading
    
    /** Empty constructor for loading */
    XMLArrayElement ()
    { 
    }
  
    /** Make a new XMLArrayElement
      * @param value This element contains an int
     */
    XMLArrayElement (int value)
    { this.value=value+"";
      type=XMLArray.INT;
    }
    
    /** Make a new XMLArrayElement
      * @param value This element contains a double
     */
    XMLArrayElement(double value)
    { this.value=value+"";
      type=XMLArray.DOUBLE;
    }
    
    /** Make a new XMLArrayElement
      * @param value This element contains a float
     */
    XMLArrayElement(float value)
    { this.value=value+"";
      type=XMLArray.FLOAT;
    }

    /** Make a new XMLArrayElement
      * @param value This element contains a boolean
     */
    XMLArrayElement(boolean value)
    { this.value=value+"";
      type=XMLArray.BOOLEAN;
    }
    
    /** Make a new XMLArrayElement
      * @param value This element contains an Integer (object)
     */
    XMLArrayElement(Integer value)
    { nullValue=(value==null);
      if (! nullValue)
         this.value=value.intValue()+"";
      type=XMLArray.INT_CLASS;
    }

    /** Make a new XMLArrayElement
      * @param value This element contains a Float (object)
     */
    XMLArrayElement(Float value)
    { nullValue=(value==null);
      if (! nullValue)
         this.value=value.floatValue()+"";
      type=XMLArray.FLOAT_CLASS;
    }

    /** Make a new XMLArrayElement
      * @param value This element contains a Double (object)
     */
    XMLArrayElement(Double value)
    { nullValue=(value==null);
      if (! nullValue)
         this.value=value.doubleValue()+"";
      type=XMLArray.DOUBLE_CLASS;
    }

    /** Make a new XMLArrayElement
      * @param value This element contains a Boolean (object)
     */
    XMLArrayElement(Boolean value)
    { nullValue=(value==null);
      if (! nullValue)
         this.value=value.booleanValue()+"";
      type=XMLArray.BOOLEAN_CLASS;
    }
    
    /** Make a new XMLArrayElement
      * @param value This element contains a String
     */
    XMLArrayElement(String value)
    { nullValue=(value==null);
      if (! nullValue)
         this.value=XMLUtils.toXMLString(value);
      type=XMLArray.STRING;
    }

    /** Make a new XMLArrayElement
      * @param value This element contains a XMLSerializable object
     */
    XMLArrayElement(XMLSerializable value)
    { nullValue=(value==null);
      if (! nullValue)
         this.value=value.getClass().getName();
      type=XMLArray.XMLSERIALIZABLE;
      child=value;
    }

    /** Make a new XMLArrayElement
      * @param value This element contains another XMLArray
     */
    XMLArrayElement(XMLArray value)
    { nullValue=(value==null);
      if (! nullValue)
         this.value=value.getClass().getName();
      type=XMLArray.OTHER_ARRAY;
      child=value;
    } 

    /** Make a new XMLArrayElement
      * @param value This element contains a byte 
     */
    XMLArrayElement(byte value)
    { this.value=(new Byte(value)).toString();
      type=XMLArray.BYTE;
    }
	 
    /** Make a new XMLArrayElement
      * @param value This element contains a Byte object
     */
    XMLArrayElement(Byte value)
    { nullValue=(value==null);
      if (! nullValue)
         this.value=value.toString();
      type=XMLArray.BYTE_CLASS;
    }
    
    // Get result methods
    
    /** @returns The int which this XMLArrayElement contains
      * @throws XMLInvalidInputException If this XMLArrayElement doesn't
      *         contain an int value
     */
    int getIntValue() throws XMLInvalidInputException
    { assertType(XMLArray.INT);
      return Integer.parseInt(value);
    }
    
    /** @returns The double which this XMLArrayElement contains
      * @throws XMLInvalidInputException If this XMLArrayElement doesn't
      *         contain a double value
     */
    double getDoubleValue() throws XMLInvalidInputException
    { assertType(XMLArray.DOUBLE);
      return Double.parseDouble(value);
    }

    /** @returns The float which this XMLArrayElement contains
      * @throws XMLInvalidInputException If this XMLArrayElement doesn't
      *         contain a float value
     */
    float getFloatValue() throws XMLInvalidInputException
    { assertType(XMLArray.FLOAT);
      return Float.parseFloat(value);
    }

    /** @returns The boolean which this XMLArrayElement contains
      * @throws XMLInvalidInputException If this XMLArrayElement doesn't
      *         contain an boolean value
     */
    boolean getBooleanValue() throws XMLInvalidInputException
    { assertType(XMLArray.BOOLEAN);
      return (new Boolean(value)).booleanValue();
    }
    
    /** @returns The Integer object which this XMLArrayElement contains
      * @throws XMLInvalidInputException If this XMLArrayElement doesn't
      *         contain an Integer object value
     */
    Integer getIntClassValue () throws XMLInvalidInputException
    { assertType(XMLArray.INT_CLASS);
      if (nullValue)
         return null;
      else 
         return new Integer(value);
    }
    
    /** @returns The Double object which this XMLArrayElement contains
      * @throws XMLInvalidInputException If this XMLArrayElement doesn't
      *         contain a Double object value
     */
    Double getDoubleClassValue () throws XMLInvalidInputException
    { assertType(XMLArray.DOUBLE_CLASS);
      if (nullValue)
         return null;
      else 
         return new Double(value);
    }

    /** @returns The Float object which this XMLArrayElement contains
      * @throws XMLInvalidInputException If this XMLArrayElement doesn't
      *         contain a Float object value
     */
    Float getFloatClassValue () throws XMLInvalidInputException
    { assertType(XMLArray.FLOAT_CLASS);
      if (nullValue)
         return null;
      else 
         return new Float(value);
    }

    /** @returns The Boolean object which this XMLArrayElement contains
      * @throws XMLInvalidInputException If this XMLArrayElement doesn't
      *         contain a Boolean object value
     */
    Boolean getBooleanClassValue () throws XMLInvalidInputException
    { assertType(XMLArray.BOOLEAN_CLASS);
      if (nullValue)
         return null;
      else 
         return new Boolean(value);
    }

    /** @returns The byte which this XMLArrayElement contains
      * @throws XMLInvalidInputException If this XMLArrayElement doesn't
      *         contain an byte value
     */
    byte getByteValue() throws XMLInvalidInputException
    { assertType(XMLArray.BYTE);
      return (new Byte(value)).byteValue();
    }
	 
    /** @returns The Byte object which this XMLArrayElement contains
      * @throws XMLInvalidInputException If this XMLArrayElement doesn't
      *         contain a Byte object value
     */
    Byte getByteClassValue () throws XMLInvalidInputException
    { assertType(XMLArray.BYTE_CLASS);
      if (nullValue)
         return null;
      else 
         return new Byte(value);
    }
    
    /** @returns The Object which this XMLArrayElement contains
      * @throws XMLInvalidInputException If this XMLArrayElement doesn't
      *         contain an Object value
     */
    Object getObjectValue () 
    { if (nullValue)
         return null;
      else 
         return result;
    }
    
    // Utility methods
    /** Make sure that the type of this XMLArrayElement matches the parameter
      * type
      * @param expectedType The type which must be equal to the the type of this
      *        XMLArrayElement
      * @throws XMLInvalidInputException If the parameter type doesn't match
      *         the type of this XMLArrayElement
     */
    void assertType (int expectedType) throws XMLInvalidInputException
    { if ( type != expectedType )
         throw new XMLInvalidInputException
	 ("Problem with extracting result from a XMLArrayElement :"+
	  "Expected element type "+expectedType+" is different from actual "+
	  "type "+type);
    }
    
    /** Method for creating instances of a certain XMLSerializable class */
    protected XMLSerializable getNewInstance (Class request) 
    throws ClassNotFoundException,InstantiationException,ClassCastException,
           XMLTreeException,IllegalAccessException // Lots of things that can go wrong
    { XMLSerializable result;
      if (watson==null || XMLArray.class.equals(request))
      { result=(XMLSerializable)request.newInstance();
      }
      else if (watson.canCreateInstance(request))
      { result=(XMLSerializable)watson.createInstance(request);
      }
      else
      { throw new ClassNotFoundException
        ("Problem in XMLArrayElement : my InstantiationAssistant can't create "+
	 "instances of "+request.getName());
      }
		if ( XMLArray.class.equals(request) && watson!=null)
			((XMLArray)(result)).setInstantiationAssistant(watson);
      result.setParentName(getXMLName());
      return result;
    }
				    
    // XMLSerializable implementation
    
	public void load (XMLElement myElement,XMLLoader loader) throws XMLTreeException,IOException,XMLInvalidInputException
	{ 	type=myElement.getAttribute("type").getIntValue();
	   value=myElement.getAttribute("value").getValue();
		nullValue=myElement.getAttribute("null-value").getBoolValue();
		if ( ! nullValue && 
		    (type==XMLArray.OTHER_ARRAY||type==XMLArray.XMLSERIALIZABLE))
		{  try
		   { Class childClass=Class.forName(value);
		     child=getNewInstance(childClass);
		   }
		   catch (Exception e)
		   { throw new XMLInvalidInputException
		     ("Something went wrong while trying to instantiate a new"+
		      " class in a XMLArrayElement :"+e);
		   }
		   loader.load(this,child);
		   if (type==XMLArray.OTHER_ARRAY)
		   	result=((XMLArray)(child)).getResult();
		   else
		   	result=child;
		}
	}

	public XMLElement saveSelf () throws XMLCannotSaveException
	{ 	XMLElement result=new XMLElement(XMLUtils.getLastName(this));
		result.addAttribute(new XMLAttribute("type",type));
		result.addAttribute(new XMLAttribute("value",value));
		result.addAttribute(new XMLAttribute("null-value",nullValue));
	  	return result;
	}
  
	public void saveChilds (XMLSaver saver) throws XMLTreeException,IOException,XMLCannotSaveException
	{ 	if (child!=null)
	           saver.saveObject(child);
	}

	public String getXMLName ()
	{ 	return parentName+".array.element";
	}
	
	public void setParentName (String parentName) throws XMLTreeException
	{	throw new XMLTreeException
		("XMLArrayElement does not support setParentName. It uses the "+
		 "parentname of its parent XMLArray.");
	}
	
  }
  
  
}
