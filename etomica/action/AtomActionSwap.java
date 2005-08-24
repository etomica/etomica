/*
 * Created on Jul 11, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.action;

import etomica.atom.AtomPair;
import etomica.atom.AtomSet;


public class AtomActionSwap extends AtomsetActionAdapter {

    public AtomActionSwap(AtomPair atomPair) {
        swappedPair = atomPair;
    }
    
    public void actionPerformed(AtomSet a) {
        swappedPair.atom0 = ((AtomPair) a).atom1;
        swappedPair.atom1 = ((AtomPair) a).atom0;
        wrappedAction.actionPerformed(swappedPair);
    }

    public AtomsetAction wrappedAction;
    private final AtomPair swappedPair;
}