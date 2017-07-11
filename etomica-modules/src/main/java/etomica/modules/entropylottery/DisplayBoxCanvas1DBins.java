/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.entropylottery;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.util.Iterator;

import etomica.action.activity.Controller;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.IData;
import etomica.data.IEtomicaDataSource;
import etomica.graphics.ColorSchemeCollective;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayCanvas;
import etomica.graphics.Drawable;
import etomica.units.Pixel;

/**
 * Class used to define canvas onto which configuration is drawn
 */
public class DisplayBoxCanvas1DBins extends DisplayCanvas {
    private double yScale = 0.4;
    private final AtomIteratorLeafAtoms atomIterator = new AtomIteratorLeafAtoms();
    private int[] atomCount;
    private IEtomicaDataSource extraDataSource;
    
    public DisplayBoxCanvas1DBins(DisplayBox _box, Controller controller) {
        super(controller);
        displayBox = _box;
        atomCount = new int[0];

        pixel = new Pixel(10);

        addComponentListener(new ComponentListener() {
            public void componentHidden(ComponentEvent e) {}
            public void componentMoved(ComponentEvent e) {}
            public void componentShown(ComponentEvent e) {}
            public void componentResized(ComponentEvent e) { refreshSize(); }});
    }
    
    protected void refreshSize() {
        Dimension dim = getSize();
        Vector boxDim = displayBox.getBox().getBoundary().getBoxSize();
        double px = (dim.width - 1)/(boxDim.getX(0)+displayBox.getPaddingSigma());
        if (pixel != null && pixel.toPixels() == px) {
            return;
        }
        setPixelUnit(new Pixel(px));
        displayBox.computeImageParameters();
    }
    
    public void setYScale(double newYScale) {
        yScale = newYScale;
    }
    
    public double getYScale() {
        return yScale;
    }
    
    /**
    * Sets the size of the display to a new value and scales the image so that
    * the box fits in the canvas in the same proportion as before.
    */
    public void scaleSetSize(int width, int height) {
        if(getBounds().width * getBounds().height != 0) {  //reset scale based on larger size change
            double ratio1 = (double)width/(double)getBounds().width;
            double ratio2 = (double)height/(double)getBounds().height;
            double factor = Math.min(ratio1, ratio2);
    //        double factor = (Math.abs(Math.log(ratio1)) > Math.abs(Math.log(ratio2))) ? ratio1 : ratio2;
            displayBox.setScale(displayBox.getScale()*factor);
            setSize(width, height);
        }
    }
            
          
    //Override superclass methods for changing size so that scale is reset with any size change  
    // this setBounds is ultimately called by all other setSize, setBounds methods
    public void setBounds(int x, int y, int width, int height) {
        if(width == 0 || height == 0) return;
        super.setBounds(x,y,width,height);
        createOffScreen(width,height);
    }
    
    /**
     * Sets a data source used to draw an extra set of "atoms".  The extra set
     * is drawn as thin blue bars.
     */
    public void setExtraData(IEtomicaDataSource newExtraDataSource) {
        extraDataSource = newExtraDataSource;
    }
    
    /**
    * doPaint is the method that handles the drawing of the box to the screen.
    * Several variables and conditions affect how the image is drawn.  First,
    * the Unit.Length.Sim class variable <code>TO_PIXELS</code> performs the conversion between simulation
    * length units (Angstroms) and pixels.  The default value is 10 pixels/Angstrom
    * reflecting the default size of the box (300 pixels by 300 pixels) and the
    * default box size (30 by 30 A).  
    * The field <code>scale</code> is a multiplicative factor that directly
    * scales up or down the size of the image; this value is adjusted automatically
    * whenever shells of periodic images are drawn, to permit the central image and all 
    * of the specified periodic images to fit in the drawing of the box.  
    *
    * @param g The graphic object to which the image of the box is drawn
    */
    public void doPaint(Graphics g) {
        if(!isVisible() || displayBox.getBox() == null) {return;}
        int w = getSize().width;
        int h = getSize().height;
            
        g.setColor(getBackground());
        g.fillRect(0,0,w,h);
        displayBox.computeImageParameters2(w, h);

        //do drawing of all drawing objects that have been added to the display
        for(Iterator iter=displayBox.getDrawables().iterator(); iter.hasNext(); ) {
            Drawable obj = (Drawable)iter.next();
            obj.draw(g, displayBox.getOrigin(), displayBox.getToPixels());
        }
            
        //Draw all atoms
        if(displayBox.getColorScheme() instanceof ColorSchemeCollective) {
            ((ColorSchemeCollective)displayBox.getColorScheme()).colorAllAtoms();
        }
        
        Box box = displayBox.getBox();
        Vector dimensions = box.getBoundary().getBoxSize();
        if (atomCount.length != (int)Math.round(dimensions.getX(0))) {
            atomCount = new int[(int)Math.round(dimensions.getX(0))];
        }
        for (int i=0; i<atomCount.length; i++) {
            atomCount[i] = 0;
        }
        atomIterator.setBox(box);
        atomIterator.reset();
        for (IAtom a = atomIterator.nextAtom(); a != null;
             a = atomIterator.nextAtom()) {
            int x = (int)Math.round(a.getPosition().getX(0)+dimensions.getX(0)*0.5-0.5);
            atomCount[x]++;
        }
        
        
        int[] origin = displayBox.getOrigin();
        g.setColor(Color.RED);
        int drawingHeight = displayBox.getDrawingHeight();
        for (int i=0; i<atomCount.length; i++) {
            int baseXP = origin[0] + (int)(displayBox.getToPixels()*i);
            int height = (int)(displayBox.getToPixels()*atomCount[i]*yScale);
            if (height > drawingHeight) {
                height = drawingHeight;
            }
            int baseYP = origin[1] + drawingHeight - height;
            g.fillRect(baseXP, baseYP, (int)displayBox.getToPixels(), height);
        }

        if (extraDataSource == null) {
            return;
        }
        
        IData extraData = extraDataSource.getData();
        if (extraData.getLength() != atomCount.length) {
            // we caught it at a bad time.  we'll call back later.
            return;
        }
        g.setColor(Color.BLUE);
        for (int i=0; i<atomCount.length; i++) {
            if (extraData.getValue(i) == 0) {
                continue;
            }
            int baseXP = origin[0] + (int)(displayBox.getToPixels()*i);
            int baseYP = origin[1] + drawingHeight - (int)(displayBox.getToPixels()*extraData.getValue(i)*yScale);
            int height = 5;
            if (baseYP + height > drawingHeight) {
                height = drawingHeight - baseYP;
            }
            else if (baseYP < origin[1]) {
                continue;
            }
            else if (baseYP + height < origin[1]) {
                height = 5-(baseYP-origin[1]);
                baseYP = origin[1];
            }
            g.fillRect(baseXP, baseYP, (int)displayBox.getToPixels(), height);
        }
    }
}
