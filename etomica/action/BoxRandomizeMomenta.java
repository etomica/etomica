package etomica.action;

import etomica.api.IAction;
import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IRandom;

/**
 * Randomizes the velocities of all the leaf atoms in an IBox based on the
 * Maxwell-Boltzmann distribution.
 * 
 * @author Andrew Schultz
 */
public class BoxRandomizeMomenta implements IAction {

    public BoxRandomizeMomenta(IBox box, IRandom random) {
        this.box = box;
        atomActionRandomizeVelocity = new AtomActionRandomizeVelocity(0, random);
    }

    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }
    
    public double getTemperature() {
        return temperature;
    }
    
    public IBox getBox() {
        return box;
    }
    
    public void actionPerformed() {
        atomActionRandomizeVelocity.setTemperature(temperature);
        IAtomSet leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            atomActionRandomizeVelocity.actionPerformed(leafList.getAtom(iLeaf));
        }
    }

    protected final AtomActionRandomizeVelocity atomActionRandomizeVelocity;
    protected final IBox box;
    protected double temperature;
}
