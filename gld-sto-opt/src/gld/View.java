
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

package gld;

import java.io.*;
import java.awt.*;
import java.util.*;
import java.awt.geom.*;
import java.awt.image.*;

import gld.infra.*;
import gld.tools.Tool;
import gld.xml.*;

/**
 *
 * A <code>View</code> is a graphical representation of an infrastructure.
 *
 * @author  Group GUI
 * @version 1.0
 */

public abstract class View extends Canvas implements Observer, XMLSerializable
{
  /** The back buffer */
  protected BufferedImage buffer;
  /** The transform used to transform coordinates from infra to buffer space */
  protected AffineTransform bufferTransform;
  /** The transform used to transform coordinates from infra to view space */
  protected AffineTransform viewTransform;
  
  /** Position of this view in the infrastructure in infra coords */
  protected Point viewPosition = new Point(0, 0);
  /** Position of the buffer in zoomed view coords */
  protected Point bufferPosition = new Point(0, 0);
  /** Position of the viewport in zoomed view coords */
  protected Point viewportPosition = new Point(0, 0);
  
  /** Size of this view in in pixels */
  protected Dimension viewSize;
  /** Size of the viewport in pixels */
  protected Dimension viewportSize;
  /** Size of the buffer in pixels */
  protected Dimension bufferSize;
  /** Size of the infrastructure */
  protected Dimension infraSize;
  
  /** Number of pixels added to viewport for buffer */
  protected int bufferOverflow = 50;
  /** Number of pixels the buffer always overflows */
  protected int bufferMinimumOverflow = 5;

  /** The current zoom factor */
  protected float zoomFactor = 1.0f;
  /** A list of all possible zoom factors */
  protected final float[] zoomFactors = { 0.25f, 0.5f, 0.75f, 1f, 1.5f, 2.0f, 2.5f };
  /** The index of the current zoom factor into the zoom factor list */
  protected int	zoomIndex = 3;
  
  /** List of all overlays */
  protected LinkedList overlays;
  
  

  public View(Dimension infras)
  {
  	setBackground(Color.white);
		overlays = new LinkedList();
		bufferTransform = new AffineTransform();
		viewTransform = new AffineTransform();
		
		viewSize = infraSize = infras;
		viewPosition = viewportPosition = bufferPosition =
			new Point(-(int)(viewSize.width / 2), -(int)(viewSize.height / 2));
		
		viewportSize = new Dimension(20, 20);
		bufferSize = new Dimension(100, 100);
		buffer = new BufferedImage(100, 100, BufferedImage.TYPE_INT_RGB);
		
		createTransform();
	}

	public View() { }












	/*============================================*/
	/* Basic GET and SET methods                  */
	/*============================================*/

	/** Returns the current zoom factor */
	public float getZoomFactor() { return zoomFactor; }
	/** Sets the current zoom factor */
	public void setZoomFactor(float z) { zoomFactor = z; }

	/** Returns the current zoom index */
	public int getZoomIndex() { return zoomIndex; }
	/** Sets the current zoom index */
	public void setZoomIndex(int i) { zoomIndex = i; }


	public Point getViewPosition() { return viewPosition; }
	public void setViewPosition(Point p) { viewPosition = p; }
	
	public Point getViewportPosition() { return viewportPosition; }
	public void setViewportPosition(Point p) { viewportPosition = p; }
	
	public Point getBufferPosition() { return bufferPosition; }
	public void setBufferPosition(Point p) { bufferPosition = p; }
	
	public Dimension getViewSize() { return viewSize; }
	public void setViewSize(Dimension d) { viewSize = d; }

	public Dimension getViewportSize() { return viewportSize; }
	public void setViewportSize(Dimension d) { viewportSize = d; }
	
	public Dimension getBufferSize() { return bufferSize; }
	public void setBufferSize(Dimension d) { bufferSize = d; }

	public Dimension getInfraSize() { return infraSize; }
	public void setInfraSize(Dimension d) { infraSize = d; }
	
	public int getBufferOverflow() { return bufferOverflow; }
	public void setBufferOverflow(int bo) { bufferOverflow = bo; }

	/**
	 * Adds an overlay to this view.
	 * An overlay has one of two types: buffer or view.
	 * A buffer overlay is painted directly on the back buffer,
	 * a view overlay on the front buffer.
	 *
	 * @param ol The overlay to add
	 */
	public void addOverlay(Overlay ol) {
		overlays.add(ol);
		if (ol.overlayType() == 1) repaint();
		else redraw();
	}
	
	/**
	 * Removes an overlay from this view.
	 *
	 * @param ol The overlay to remove.
	 */
	public void remOverlay(Overlay ol) {
		if (ol == null) return;
		overlays.remove(ol);
		if (ol.overlayType() == 1) repaint();
		else redraw();
	}










	/*============================================*/
	/* Misc                                       */
	/*============================================*/


	/** Increases the zoom factor */
	public void zoomIn() {
		if (zoomIndex < zoomFactors.length - 1) zoomTo(zoomIndex + 1);
	}
	/** Decreases the zoom factor. */
	public void zoomOut() {
		if (zoomIndex > 0) zoomTo(zoomIndex - 1);
	}

	/**
	 * Sets the zoom factor for this <code>View</code>.
	 *
	 * @param index The index of the new zoom factor.
	 */
	public void zoomTo(int index)
	{
		zoomIndex = index;
		zoomFactor = zoomFactors[zoomIndex];
		viewSize = new Dimension((int)(infraSize.width * zoomFactor), (int)(infraSize.height * zoomFactor));
		setSize(viewSize);
		createTransform();
		redraw();
	}

	/**
	 * Scrolls the viewport to the specified view point
	 *
	 * @param p The point in view coordinate space to scroll to
	 */
	public void scrollViewport(Point p) {
		setViewportPosition(p);

		Rectangle buf = new Rectangle(bufferPosition, bufferSize);
		Rectangle vp = new Rectangle(viewportPosition, viewportSize);
		vp.grow(bufferMinimumOverflow, bufferMinimumOverflow);
		
		if (!buf.contains(vp)) repositionBuffer();
		else repaint();
	}
	
	/** Resizes the viewport to the given size */
	public void resizeViewport(Dimension d) {
		viewportSize = d;
		resizeBuffer();
	}
	
	/** Resizes the infrastructure */
	public void resizeInfra(Dimension d) {
		infraSize = d;
		viewSize = new Dimension((int)(d.width * zoomFactor), (int)(d.height * zoomFactor));
		viewPosition = new Point(-(int)(d.width / 2), -(int)(d.height / 2));
		setSize(viewSize);
		createTransform();
		redraw();
	}





















	/*============================================*/
	/* Painting                                   */
	/*============================================*/

	public void update(Observable o, Object arg) {
		if (o instanceof Overlay) {
			if (((Overlay)o).overlayType() == 1) repaint();
			else redraw();
		}
		else redraw();
	}

	protected void redraw() {
		if (!isVisible()) return;
		drawBuffer();
		repaint();
	}

	private void drawBuffer()
	{
		Graphics2D g = buffer.createGraphics();
		g.setColor(Color.white);
		g.fillRect(0, 0, bufferSize.width, bufferSize.height);
		g.setTransform(bufferTransform);
		g.setColor(Color.black);
		fillBuffer(g);
		drawOverlays(g, 2);
	}

	protected abstract void fillBuffer(Graphics2D g);
	
	private void drawOverlays(Graphics2D g, int type) {
		Overlay ol = null;
		for (ListIterator it = overlays.listIterator(); it.hasNext(); ) {
			ol = (Overlay)it.next();
			if (ol.overlayType() == type) {
				try { ol.paint(g); }
				catch (GLDException e) { Controller.reportError(e); }
			}
		}
	}

	public void update(Graphics g)
	{
		paint(g);
	}

	public void paint(Graphics gr)
	{
		if (!isVisible()) return;
		Graphics2D g = (Graphics2D)gr;
		while (!g.drawImage(buffer, bufferPosition.x, bufferPosition.y, this));
		g.setTransform(viewTransform);
		drawOverlays(g, 1);
	}

	/** Repositions the internal buffer */
	private void repositionBuffer() {
		bufferPosition.x = viewportPosition.x - bufferOverflow;
		bufferPosition.y = viewportPosition.y - bufferOverflow;
		createTransform();
		redraw();
	}
	
	/** Resizes the internal buffer */
	private void resizeBuffer() {
		bufferSize.width = viewportSize.width + bufferOverflow * 2;
		bufferSize.height = viewportSize.height + bufferOverflow * 2;
		buffer = new BufferedImage(bufferSize.width, bufferSize.height, BufferedImage.TYPE_INT_RGB);
		repositionBuffer();
	}
	
	/** Returns the preferred size of this view. Used by ScrollPane. */
	public Dimension getPreferredSize() { return viewSize; }
	
	private void createTransform() {
		bufferTransform.setToTranslation(-bufferPosition.x, -bufferPosition.y);
		bufferTransform.scale(zoomFactor, zoomFactor);
		bufferTransform.translate(-viewPosition.x, -viewPosition.y);
		
		viewTransform.setToScale(zoomFactor, zoomFactor);
		viewTransform.translate(-viewPosition.x, -viewPosition.y);
	}
















	/*============================================*/
	/* LOAD AND SAVE                              */
	/*============================================*/

	public void load (XMLElement myElement,XMLLoader loader) throws XMLTreeException, IOException, XMLInvalidInputException
	{
 		int height, width;
		height = myElement.getAttribute("height").getIntValue();
		width = myElement.getAttribute("width").getIntValue();
		zoomIndex = myElement.getAttribute("zoom-index").getIntValue();
		zoomFactor = zoomFactors[zoomIndex];
		viewSize = new Dimension(width, height);
		resizeBuffer();
 	}

	public XMLElement saveSelf () throws XMLCannotSaveException
	{
		XMLElement result=new XMLElement(getXMLName());
		result.addAttribute(new XMLAttribute ("width",buffer.getWidth()));
		result.addAttribute(new XMLAttribute ("height",buffer.getHeight()));
		result.addAttribute(new XMLAttribute ("zoom-index",zoomIndex));
		return result;
	}

 	public void saveChilds (XMLSaver saver) throws XMLTreeException, IOException, XMLCannotSaveException
 	{ 	// The view doesn't have any child objects
 	}

 	public String getXMLName ()
 	{ 	return "view";
 	}
	
	public void setParentName (String parentName) throws XMLTreeException
	{	throw new XMLTreeException
	        ("Cannot set parent for XML root class view");
	}



















	/*============================================*/
	/* COORDINATE CONVERSIONS                     */
	/*============================================*/

	/** Converts the given point from infrastructure to view coordinate space */
	public Point toView(Point p) {
		return new Point(
			(int)((p.x - viewPosition.x) * zoomFactor),
			(int)((p.y - viewPosition.y) * zoomFactor)
		);
	}
	
	/** Converts the given point from view to infrastructure coordinate space */
	public Point toInfra(Point p) {
		return new Point(
			(int)(p.x / zoomFactor + viewPosition.x),
			(int)(p.y / zoomFactor + viewPosition.y)
		);
	}
}