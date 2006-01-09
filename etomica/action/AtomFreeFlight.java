/*
 * History
 * Created on Nov 18, 2004 by kofke
 */
package etomica.action;

import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.space.ICoordinateKinetic;


/**
 * Translates the atom by the amount it would move in free (ballistic) flight for a specified time interval.
 * Uses the atom's current momentum to determine this displacement.
 */
public class AtomFreeFlight extends AtomActionAdapter {
    private double tStep = 0.0;
    public void actionPerformed(Atom a) {
        ((AtomLeaf)a).coord.position().PEa1Tv1(tStep,((ICoordinateKinetic)((AtomLeaf)a).coord).velocity());
    }
    public void actionPerformed(Atom a, double t) {
        tStep = t;
        actionPerformed(a);
    }
    public void setTStep(double t) {tStep = t;}
    public double getTStep() {return tStep;}
}