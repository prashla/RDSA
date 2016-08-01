
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

import java.io.*;
/** This interface has to be implemented by objects that want to 
  * be loaded and saved via the parser.
 */
 
public interface XMLSerializable
{ 
 /** Load this XMLSerializable
  *  @param myself The XMLElement which represents this object in the
  *         XML tree. It can contain attributes which hold information
  *         about the object.
  *  @param loader The XMLLoader which this XMLSerializable can use to
  *         load child objects.
  *  @throws XMLTreeException The parser can throw this exception if it
  *          is called by the XMLSerializable. The XMLSerializable should NOT
  *          throw this exception by itself and preferrably not catch it too.
  *  @throws IOException Thrown in case of an read error in the XML file.
  *  @throws XMLInvalidInputException The XMLSerializable can throw this
  *          exception if it cannot load itself or one of its child objects for
  *          whatever reason.
  */
 void load (XMLElement myself,XMLLoader loader) throws XMLTreeException,IOException,XMLInvalidInputException;

 /** @return The XMLElement which represents the internal datastate. This
  *          function is not meant to save child objects. That is done in
  *          void saveChilds()
  *  @throws XMLCannotSaveException The XMLSerializable can throw this
  *          exception if it cannot save itself or one of its child objects for
  *          whatever reason.
  */
 XMLElement saveSelf () throws XMLCannotSaveException;
  
 /** This method gives an XMLSerializable the opportunity to save 
  *  its child objects.
  *  @param saver The XMLSaver that the XMLSerializable can use to save its
  *         child objects.
  *  @throws XMLTreeException The parser can throw this exception if it
  *          is called by the XMLSerializable. The XMLSerializable should NOT
  *          throw this exception by itself and preferrably also not catch it too.
  *  @throws IOException Thrown in case the parser cannot write to the file
  *  @throws XMLCannotSaveException The XMLSerializable can throw this
  *          exception if it cannot save itself or one of its child objects for
  *          whatever reason.
  */
 void saveChilds (XMLSaver saver) throws XMLTreeException,IOException,XMLCannotSaveException;

 /** @return The tagname of the XML element in which information about
  *         this object is stored prepended by a dot separated list 
  *         of the tagnames of the parents of this object. (e.g. 
  *         "model.infrastructure.node")
  */
 String getXMLName ();
 
 /** Sets a new parent name for this object
   * @param parentName The new parent name
   * @throws XMLTreeException If this object does not support setting other
   *         parent names.
  */
 void setParentName (String parentName) throws XMLTreeException;
}
