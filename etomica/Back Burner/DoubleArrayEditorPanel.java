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
import java.awt.event.*;
import java.text.*;
import java.lang.reflect.*;
import java.beans.*;
import java.io.*;

/**
* This class displays the values of an array of doubles
* in fields which can be edited to change specific values.
*/
public class DoubleArrayEditorPanel extends Panel implements ItemListener, Serializable {

    private PropertyEditorSupport editor;
    private double[] theArray;
    private int currentIndex = 0;
    private NumberFormat fmt = NumberFormat.getNumberInstance();
    private TextField sizeField = new TextField(4);
    private TextField valueField = new TextField(12);
    private List elementList = new List();

    /**
    * This class is instantiated by the property editor which
    * is passed as type PropertyEditorSupport
    */
    public DoubleArrayEditorPanel(PropertyEditorSupport ed) {
        editor = ed;

        if(ed == null) System.out.println("editor is send a null to constructor");
        else System.out.println("editor is not null");

        if(ed.getValue() == null) System.out.println("ed.getValue is returning a null");
        else System.out.println("ed.getValue is not null");

        double[] myArray = new double[2];
        myArray = (double[])editor.getValue();

        setArray((double[])ed.getValue());
        System.out.println("another hello from DoubleArrayEditorPanel Constructor, length = " + theArray.length);
        setLayout(new GridBagLayout());
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.weightx = 0;
        gbc.weighty = 0;
        gbc.fill = GridBagConstraints.NONE;
        gbc.anchor = GridBagConstraints.EAST;
        add(new Label("Size"), gbc, 0, 0, 1, 1);
        add(new Label("Elements"), gbc, 0, 1, 1, 1);
        gbc.weightx = 100;
        gbc.anchor = GridBagConstraints.WEST;
        add(sizeField, gbc, 1, 0, 1, 1);
        gbc.fill = GridBagConstraints.HORIZONTAL;
        add(valueField, gbc, 1, 1, 1, 1);
        gbc.weighty = 100;
        gbc.fill = GridBagConstraints.BOTH;
        add(elementList, gbc, 1, 2, 1, 1);
        gbc.fill = GridBagConstraints.NONE;

        /**
        * Adds this object as an ItemListener for the property
        * elementList
        */
        elementList.addItemListener(this);

        /**
        * Sets the length of the array.
        * Overrides method addKeyListener of sizeField which holds length
        * of array.  If enter key is pressed, calls resizeArray method.
        * @parm anonymous class of type KeyAdapter of which method
        * keyPressed is overridden to call resizeArray if the enter key
        * is pressed.
        * @see #resizeArray()
        */
        sizeField.addKeyListener(new KeyAdapter() {
            public void keyPressed(KeyEvent evt) {
                if (evt.getKeyCode() == KeyEvent.VK_ENTER) {
                    resizeArray();
                }
            }
            }
        );

        /**
        * Sets the length of the array.
        * Overrides method addFocusListener of sizeField which holds length
        * of array.  If lost focus is not temporary, calls resizeArray method.
        * @parm anonymous class of type KeyAdapter of which method
        * focusLost is overridden to call resizeArray if the lost focus is not
        * temporary, (i.e. a popup menu).
        * @see #resizeArray()
        */
        sizeField.addFocusListener(new FocusAdapter() {
            public void focusLost(FocusEvent evt) {
                if (!evt.isTemporary()) {
                    resizeArray();
                }
            }
        });

        /**
        * Sets the value of an element in the array.
        * Overrides method addKeyListener of valueField which holds a particular
        * value in the array.  If enter key is pressed, calls changeValue method.
        * @parm anonymous class of type KeyAdapter of which method
        * keyPressed is overridden to call changeValue if the enter key
        * is pressed.
        * @see #changeValue()
        */
        valueField.addKeyListener(new KeyAdapter() {
            public void keyPressed(KeyEvent evt) {
                if (evt.getKeyCode() == KeyEvent.VK_ENTER) {
                    changeValue();
                }
            }
        });

        /**
        * Sets the value of an element in the array.
        * Overrides method addFocusListener of valueField which holds a
        * particular value in the array.  If focus is lost, calls changeValue method.
        * @parm anonymous class of type FocusAdapter of which method
        * focusLost is overridden to call changeValue if lost focus is not
        * temporary (i.e. popup menu is a tempory loss of focus).
        * @see #changeValue()
        */
        valueField.addFocusListener(new FocusAdapter() {
            public void focusLost(FocusEvent evt) {
                if (!evt.isTemporary()) {
                    changeValue();
                }
            }
        });
    }

    /**
    * Used to add components elementList, sizeField and valueField
    */
    public void add(Component c, GridBagConstraints gbc,
      int x, int y, int w, int h) {
        gbc.gridx = x;
        gbc.gridy = y;
        gbc.gridwidth = w;
        gbc.gridheight = h;
        add(c, gbc);
    }

    /**
    * Resizes the array of doubles by calling arrayGrow
    * and sets the array of the property editor.
    * @see #arrayGrow
    */
    public void resizeArray() {
        fmt.setParseIntegerOnly(true);
        int s = 0;
        try {
            s = fmt.parse(sizeField.getText()).intValue();
         if (s < 0)
            throw new ParseException("Out of bounds", 0);
        }
        catch(ParseException e) {
            sizeField.requestFocus();
            return;
        }
        if (s == theArray.length) return;
        setArray((double[])arrayGrow(theArray, s));
        editor.setValue(theArray);
        editor.firePropertyChange();
        System.out.println("Hello from DoubleArrayEditorPanel.resize");
    }

    /**
    * Sets the value of a specific element of the array and
    * changes the corresponding value in the propery editor.
    */
    public void changeValue() {
        double v = 0;
        fmt.setParseIntegerOnly(false);
        try {
            v = fmt.parse(valueField.getText()).doubleValue();
        }
        catch(ParseException e) {
            valueField.requestFocus();
            return;
        }
        setArray(currentIndex, v);
        editor.firePropertyChange();
    }

    /**
    * This method is called by elementlist which is a class of
    * type List.  List implements ItemSelectable which fires
    * processItemEvent (and so calling this method) anytime one
    * of its "items" has changed
    * @parm evt ItemEvent
    */
    public void itemStateChanged(ItemEvent evt) {
        if (evt.getStateChange() == ItemEvent.SELECTED) {
            int i = elementList.getSelectedIndex();
            valueField.setText("" + theArray[i]);
            currentIndex = i;
        }
    }

    /**
    * Changes the length of the array.
    */
    
    static Object arrayGrow(Object a, int newLength) {
        Class cl = a.getClass();
        if (!cl.isArray()) return null;
        Class componentType = a.getClass().getComponentType();
        int length = Array.getLength(a);

        Object newArray = Array.newInstance(componentType,
         newLength);
        System.arraycopy(a, 0, newArray, 0,
         Math.min(length, newLength));
        return newArray;
    }

    /**
    * @return array as an array of doubles
    */
    public double[] getArray() {
        return (double[])theArray.clone();
    }

    /**
    * @parm v sets values of array to that of v.
    */
    public void setArray(double[] v) {
        if (v == null) theArray = new double[0];
        else theArray = v;
        sizeField.setText("" + theArray.length);
        elementList.removeAll();
        for (int i = 0; i < theArray.length; i++)
            elementList.add("[" + i + "] " + theArray[i]);
        if (theArray.length > 0) {
            valueField.setText("" + theArray[0]);
            elementList.select(0);
            currentIndex = 0;
        }
        else
            valueField.setText("");
    }

    /**
    * @return array[i]
    */
    public double getArray(int i) {
        if (0 <= i && i < theArray.length) return theArray[i];
        return 0;
    }

    /**
    * @parm i, index of element to set.
    * @parm value, value of ith array element to be set.
    */
    public void setArray(int i, double value) {
        if (0 <= i && i < theArray.length) {
            theArray[i] = value;
            elementList.replaceItem("[" + i + "] " + value, i);
            int previous = elementList.getSelectedIndex();
            elementList.select(i);
            valueField.setText("" + value);
            elementList.select(previous);
        }
    }
}
