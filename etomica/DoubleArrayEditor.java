/*
 * Cay S. Horstmann & Gary Cornell, Core Java
 * Published By Sun Microsystems Press/Prentice-Hall
 * Copyright (C) 1997 Sun Microsystems Inc.
 * All Rights Reserved.
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for NON-COMMERCIAL purposes
 * and without fee is hereby granted provided that this
 * copyright notice appears in all copies.
 *
 * THE AUTHORS AND PUBLISHER MAKE NO REPRESENTATIONS OR
 * WARRANTIES ABOUT THE SUITABILITY OF THE SOFTWARE, EITHER
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 * PARTICULAR PURPOSE, OR NON-INFRINGEMENT. THE AUTHORS
 * AND PUBLISHER SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED
 * BY LICENSEE AS A RESULT OF USING, MODIFYING OR DISTRIBUTING
 * THIS SOFTWARE OR ITS DERIVATIVES.
 */

/**
 * @version 1.10 27 Oct 1997
 * @author Cay Horstmann
 */
package simulate;
import java.awt.*;
import java.beans.*;

public class DoubleArrayEditor extends PropertyEditorSupport implements java.io.Serializable {
    /**
    * Saves the values entered at design time so they are consistent
    * between edits.
    */
    static double[] array;
    /**
    * Called by the Bean Builder to instantiate DoubleArrayEditorPanel
    * @parm this send this instance of the Double Array Editor
    * @return a DoubleArrayEditorPanel to the Bean Builder
    * and is used to allow design-time display and edit
    */
    public Component getCustomEditor() {
        return new DoubleArrayEditorPanel(this);
    }

    /**
    * overrides method of the interface PropertyEditor and
    * indicates there is a custom editor (DoubleArrayEditorPanel)
    * @return true
    */
    public boolean supportsCustomEditor() {
        return true;
    }

    /**
    * overrides method of the interface PropertyEditor and
    * indicates this property editor does not support "getting"
    * and "setting" this property as a text string.
    * @return null
    */
    public String getAsText(){
        return null;
    }

    /**
    * overrides method of the interface PropertyEditor and
    * indicates this property editor is able to "paint" the
    * property.
    * @return true
    */
    public boolean isPaintable(){
        return true;
    }

    /**
    * overrides method of the interface PropertyEditor and
    * Implementation of how the property is "painted"
    * @parm Graphics
    * @parm Rectangle
    */
    public void paintValue(Graphics g, Rectangle box) {
        double[] values = (double[]) getValue();
        String s = "";
        for (int i = 0; i < 3; i++){
            if (values.length > i) s = s + values[i];
            if (values.length > i + 1) s = s + ", ";
        }
        if (values.length > 3) s += "...";
        g.setColor(Color.white);
        g.fillRect(box.x, box.y, box.width, box.height);
        g.setColor(Color.black);
        /*
        FontMetrics is a convenience class used to get various
        dimesions which are particular to specific fonts.  Used
        to display property within given space
        */
        FontMetrics fm = g.getFontMetrics();
        int w = fm.stringWidth(s);
        int x = box.x;
        if (w < box.width) x += (box.width - w) / 2;
        int y = box.y + (box.height - fm.getHeight()) / 2 + fm.getAscent();
        g.drawString(s, x, y);
    }

    /**
    * Used by Bean Builder to write a string to the applet/application which
    * initializes the properties at run time to the values chosen at Design Time
    * @return String
    */
    public String getJavaInitializationString() {
        //Bean Builder uses this to write string to applet.
        //array = (double[])getValue();
        String s = "new double[]{" + array[0];
 //       System.out.println("From InitializaionString " + s);
        for(int i = 1; i < array.length; i++){
            s += ", " + array[i];
        }
        s += "}";
 //       System.out.println("from InitializaionString, after for loop " + s);
        return s;
    }

    /**
    * Sets the new value of the property "array"
    * @parm value
    */
    public void setValue(Object value) {
        array = (double[])value;
    }

    /**
    * @return array
    */
    public Object getValue(){
        return array;
    }
/*
    public synchronized void addPropertyChangeListener(PropertyChangeListener listener) {
        System.out.println("Hello from DoubleArrayEditor.addPropertyChangeListener");
        addPropertyChangeListener(listener);
    }
    */
}


