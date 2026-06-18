/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import java.awt.*;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;

/**
 * Superclass of all classes that display something from the simulation.
 * Included are displays of graphical and tabular data, and views of the
 * molecules as they move about during the simulation.
 *
 * @author David Kofke
 */
public abstract class Display {

    public Display() {
    }

    /**
     * Method of Simulation.GraphicElement interface. Default action is to
     * return this Display (which is a Panel) as the graphic object. May
     * override in subclass to return a more appropriate graphical element, or
     * none at all.
     */
    public abstract Component graphic();

    /**
     * Accessor method of the label describing the display.
     * 
     * @return The given label
     */
    public String getLabel() {
        return label;
    }

    /**
     * Accessor method of the label describing the display.
     * 
     * @param label
     *            The label string describing the display
     */
    public void setLabel(String label) {
        String oldLabel = this.label;
        this.label = label;
        support.firePropertyChange("label", oldLabel, label);
    }

    /**
     * Overrides the Object class toString method to have it return the output
     * of getName
     */
    public String toString() {
        if (label == null) {
            return label+" "+this.getClass().getName();
        }
        return this.getClass().getName();
    }

    public void addPropertyChangeListener(PropertyChangeListener listener) {
        support.addPropertyChangeListener(listener);
    }

    public void removePropertyChangeListener(PropertyChangeListener listener) {
        support.removePropertyChangeListener(listener);
    }

    protected PropertyChangeSupport support = new PropertyChangeSupport(this);
    protected String label;
}