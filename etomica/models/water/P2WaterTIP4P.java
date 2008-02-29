
package etomica.models.water;

import etomica.api.IVector;
import etomica.atom.AtomSet;
import etomica.atom.IAtomPositioned;
import etomica.atom.IMolecule;
import etomica.box.Box;
import etomica.potential.Potential2;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.units.Electron;
import etomica.units.Kelvin;

/** 
 * 
 * @author kofke
 *
 * TIP4P potential for water.  Requires the molecule node to be an
 * AtomTreeNodeWaterTIP4P.  This does not apply periodic boundary conditions.
 */
public class P2WaterTIP4P extends Potential2 {

    public P2WaterTIP4P(Space space) {
	    super(space);
	    work = space.makeVector();
	    shift = space.makeVector();
        sigma = 3.1540;
        sigma2 = sigma*sigma;
        epsilon = Kelvin.UNIT.toSim(78.02);
        epsilon4 = 4*epsilon;
        chargeH = Electron.UNIT.toSim(0.52);
        chargeM = Electron.UNIT.toSim(-1.04);
        chargeHH = chargeH*chargeH;
        chargeHM = chargeM*chargeH;
        chargeMM = chargeM*chargeM;
    }   

    public double energy(AtomSet atoms){
        double sum = 0.0;
        double r2 = 0.0;

        AtomSet water1Atoms = ((IMolecule)atoms.getAtom(0)).getChildList();
        AtomSet water2Atoms = ((IMolecule)atoms.getAtom(1)).getChildList();

        IVector O1r = ((IAtomPositioned)water1Atoms.getAtom(SpeciesWater4P.indexO)).getPosition();
        IVector O2r = ((IAtomPositioned)water2Atoms.getAtom(SpeciesWater4P.indexO)).getPosition();
        
        work.Ev1Mv2(O1r, O2r);
        shift.Ea1Tv1(-1,work);
        boundary.nearestImage(work);
        shift.PE(work);
        final boolean zeroShift = shift.squared() < 0.1;
        r2 = work.squared();

        if(r2<1.2) return Double.POSITIVE_INFINITY;

        double s2 = sigma2/(r2);
        double s6 = s2*s2*s2;
        sum += epsilon4*s6*(s6 - 1.0);
        
        IVector H11r = ((IAtomPositioned)water1Atoms.getAtom(SpeciesWater4P.indexH1)).getPosition();
        IVector H12r = ((IAtomPositioned)water1Atoms.getAtom(SpeciesWater4P.indexH2)).getPosition();
        IVector H21r = ((IAtomPositioned)water2Atoms.getAtom(SpeciesWater4P.indexH1)).getPosition();
        IVector H22r = ((IAtomPositioned)water2Atoms.getAtom(SpeciesWater4P.indexH2)).getPosition();

        IVector M1r = ((IAtomPositioned)water1Atoms.getAtom(SpeciesWater4P.indexM)).getPosition();
        IVector M2r = ((IAtomPositioned)water1Atoms.getAtom(SpeciesWater4P.indexM)).getPosition();

        if (zeroShift) {
            r2 = M1r.Mv1Squared(H21r);
            sum += chargeHM/Math.sqrt(r2);
            r2 = M1r.Mv1Squared(H22r);
            sum += chargeHM/Math.sqrt(r2);
            r2 = H11r.Mv1Squared(M2r);
            sum += chargeHM/Math.sqrt(r2);
            r2 = H11r.Mv1Squared(H21r);
            sum += chargeHH/Math.sqrt(r2);
            r2 = H11r.Mv1Squared(H22r);
            sum += chargeHH/Math.sqrt(r2);
            r2 = H12r.Mv1Squared(M2r);
            sum += chargeHM/Math.sqrt(r2);
            r2 = H12r.Mv1Squared(H21r);
            sum += chargeHH/Math.sqrt(r2);
            r2 = H12r.Mv1Squared(H22r);
            sum += chargeHH/Math.sqrt(r2);
        }
        else {
            shift.PE(M1r);
            r2 = H21r.Mv1Squared(shift);
            shift.ME(M1r);
            sum += chargeHM/Math.sqrt(r2);

            shift.PE(M1r);
            r2 = H22r.Mv1Squared(shift);
            shift.ME(M1r);
            sum += chargeHM/Math.sqrt(r2);

            shift.PE(H11r);
            r2 = M2r.Mv1Squared(shift);
            shift.ME(H11r);
            sum += chargeHM/Math.sqrt(r2);

            shift.PE(H11r);
            r2 = H21r.Mv1Squared(shift);
            shift.ME(H11r);
            sum += chargeHH/Math.sqrt(r2);

            shift.PE(H11r);
            r2 = H22r.Mv1Squared(shift);
            shift.ME(H11r);
            sum += chargeHH/Math.sqrt(r2);

            shift.PE(H12r);
            r2 = M2r.Mv1Squared(shift);
            shift.ME(H12r);
            sum += chargeHM/Math.sqrt(r2);

            shift.PE(H12r);
            r2 = H21r.Mv1Squared(shift);
            shift.ME(H12r);
            sum += chargeHH/Math.sqrt(r2);

            shift.PE(H12r);
            r2 = H22r.Mv1Squared(shift);
            shift.ME(H12r);
            sum += chargeHH/Math.sqrt(r2);
        }

        return sum;																					        
    }

    public double getSigma() {return sigma;}

    public final double getRange() {
        return Double.POSITIVE_INFINITY;
    }
    
    public double getEpsilon() {return epsilon;}
    
    public void setBox(Box box) {
        boundary = box.getBoundary();
    }

    private static final long serialVersionUID = 1L;
    protected final double sigma , sigma2;
    protected final double epsilon, epsilon4;
    protected Boundary boundary;
    protected final IVector work, shift;
    protected final double chargeH;
    protected final double chargeM;
    protected final double chargeHH;
    protected final double chargeHM;
    protected final double chargeMM;
}