package etomica.action;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;


/**
 * Translates the atom by the amount it would move in free (ballistic) flight for a specified time interval.
 * Uses the atom's current momentum to determine this displacement.
 */
public class AtomFreeFlight implements AtomAction {
    private static final long serialVersionUID = 1L;
    private double tStep = 0.0;
    public void actionPerformed(IAtom a) {
        IAtomKinetic atom = (IAtomKinetic)a;
        atom.getPosition().PEa1Tv1(tStep,atom.getVelocity());
    }
    public void actionPerformed(IAtom a, double t) {
        tStep = t;
        actionPerformed(a);
    }
    public void setTStep(double t) {tStep = t;}
    public double getTStep() {return tStep;}
}