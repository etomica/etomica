/*
 * History
 * Created on Nov 18, 2004 by kofke
 */
package etomica.action;

import etomica.Atom;


/**
 * Translates the atom by the amount it would move in free (ballistic) flight for a specified time interval.
 * Uses the atom's current momentum to determine this displacement.
 */
public class AtomFreeFlight extends AtomActionAdapter {
    private double tStep = 0.0;
    public void actionPerformed(Atom a) {
        if(a.coord.isStationary()) {return;}  //skip if atom is stationary
        a.coord.freeFlight(tStep);  // r += tStep*p/m
    }
    public void actionPerformed(Atom a, double t) {
        tStep = t;
        actionPerformed(a);
    }
    public void setTStep(double t) {tStep = t;}
    public double getTStep() {return tStep;}
}