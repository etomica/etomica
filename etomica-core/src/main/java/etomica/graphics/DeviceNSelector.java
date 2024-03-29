/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import etomica.action.ActionGroupSeries;
import etomica.action.IAction;
import etomica.action.controller.Controller;
import etomica.box.Box;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierNMolecule;
import etomica.species.ISpecies;

/**
 * Slider that selects the number of atoms of a given species in a box.
 *
 * @author David Kofke
 */
public class DeviceNSelector extends DeviceSlider {
    
    public DeviceNSelector() {
        this(null);
    }
    
    public DeviceNSelector(Controller controller) {
        super(controller);
    }

    public void setEnabled(boolean newIsEnabled) {
        slider.setEnabled(newIsEnabled);
        textField.setEnabled(newIsEnabled);
    }

    public void setResetAction(IAction newResetAction) {
        resetAction = newResetAction;
        if (modifyAction != null) {
            targetAction = new ActionGroupSeries(modifyAction, resetAction);
        }
    }

    /**
     * Returns the action used to "reset" the simulation after changing the 
     * number of molecules, SimulationRestart by default.
     */
    public IAction getResetAction() {
        return resetAction;
    }

    public void setBox(Box newBox) {
        box = newBox;
        if (species != null) {
            init();
        }
    }
    
    public void setSpecies(ISpecies newSpecies) {
        species = newSpecies;
        if (box != null) {
            init();
        }
    }
    
    public Box getBox() {
        return box;
    }
    
    public ISpecies getSpecies() {
        return species;
    }
    
    public void setModifier(Modifier newModifier) {
        super.setModifier(newModifier);
        if (resetAction != null) {
            targetAction = new ActionGroupSeries(modifyAction, resetAction);
        }
    }
    
    protected void init() {
        setMinimum(0);
        int max = 60;
        int nMolecules = box.getNMolecules(species);
        if (nMolecules > max) max = nMolecules;
        setMaximum(max);
        slider.setSnapToTicks(false);
        slider.setMajorTickSpacing(10);
        graphic().setSize(new java.awt.Dimension(40, 30));
        setModifier(new ModifierNMolecule(box, species));

        setLabel("Number of molecules");
    }
    
    protected IAction resetAction;
    protected ISpecies species;
    protected Box box;
} //end of DeviceNSelector
