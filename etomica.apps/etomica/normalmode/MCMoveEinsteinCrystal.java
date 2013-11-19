package etomica.normalmode;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IRandom;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.iterator.AtomIterator;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.space.ISpace;

/**
 * MC move whose purpose in life is to sample an  Einstein crystal.
 * Since the energy is harmonic, each configuration can be
 * independent.
 *
 * @author Andrew Schultz
 */
public class MCMoveEinsteinCrystal extends MCMoveBox {

    public MCMoveEinsteinCrystal(ISpace space, IRandom random) {
        super(null);
        this.random = random;
        fixedCOM = true;
        dr = space.makeVector();
    }
    
    public void setFixedCOM(boolean newFixedCOM) {
        this.fixedCOM = newFixedCOM;
    }
    
    public boolean getFixedCOM() {
        return fixedCOM;
    }

    /**
     * Sets the Einstein spring constant
     */
    public void setAlphaEin(double newAlpha) {
        alpha = newAlpha;
    }
    
    public double getAlpha() {
        return alpha;
    }
    
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }
    
    public void setCoordinateDefinition(CoordinateDefinition newCoordinateDefinition) {
        coordinateDefinition = newCoordinateDefinition;
    }
    
    public boolean doTrial() {
        IAtomList atomList = box.getLeafList();
        double einFac = Math.sqrt(temperature/alpha);
        int end = atomList.getAtomCount();
        if (fixedCOM) {
            // each motion (i) also acts on all particles but in the opposite direction and less by 1/N (so that COM is unchanged) 
            // for normalization, each motion is a bit more by 1/sqrt(1-1/N) = sqrt(N/(N-1))
            einFac = Math.sqrt(temperature/alpha/Math.sqrt(1.0-1.0/end));
            // don't expclitly move the last atom
            end--;
            dr.E(0);
        }
        for (int i=0; i<end; i++) {
            IAtom a = atomList.getAtom(i);
            IVectorMutable p = a.getPosition();
            IVector site = coordinateDefinition.getLatticePosition(a);
            for (int k=0; k<p.getD(); k++) {
                p.setX(k, einFac * random.nextGaussian() );
            }
            if (fixedCOM) {
                dr.ME(p);
            }
            p.PE(site);
        }
        if (fixedCOM) {
            dr.TE(1.0/(end+1));
            for (int i=0; i<=end; i++) {
                IAtom a = atomList.getAtom(i);
                IVectorMutable p = a.getPosition();
                p.ME(dr);
            }
        }
        return true;
    }

    public AtomIterator affectedAtoms() {return null;}

    public double energyChange() {return 0;}

    public double getA() {return 1;}

    public double getB() {return 0;}

    public void acceptNotify() {}

    public void rejectNotify() {}

    protected double alpha;
    protected double temperature;
    protected final IRandom random;
    protected CoordinateDefinition coordinateDefinition;
    protected boolean fixedCOM;
    protected final IVectorMutable dr;
}
