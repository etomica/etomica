package etomica.modules.rheology;

import etomica.api.IAtomList;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.IVectorMutable;
import etomica.integrator.IntegratorMD;
import etomica.space.ISpace;

/**
 * Integrator for Brownian dynamics of a polymer in a flow field.
 *
 * @author Andrew Schultz
 */
public class IntegratorPolymer extends IntegratorMD {

    public IntegratorPolymer(IPotentialMaster potentialMaster, IRandom random,
            double timeStep, double temperature, ISpace _space) {
        super(potentialMaster, random, timeStep, temperature, _space);
        center = _space.makeVector();
        drPrev = _space.makeVector();
        drNext = _space.makeVector();
    }

    /* no point in trying to thermostat */
    public void doThermostat() {}

    public void doStepInternal() {
        super.doStepInternal();

        double z = 1;
        if (a > 0) {
            z = (1-a*a)/(1+a);
        }
        double sqa = a>0 ? Math.sqrt(a) : 0;
        double srdt = shearRate*timeStep;
        if (sqa > 1.0/160.0) {
            srdt /= 162*sqa;
        }
        IMoleculeList molecules = box.getMoleculeList();
        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            center.E(0);
            IMolecule molecule = molecules.getMolecule(i);
            IAtomList atoms = molecule.getChildList();
            IVectorMutable p0 = atoms.getAtom(0).getPosition();
            IVectorMutable p1 = atoms.getAtom(1).getPosition();
            //calculate 0-1 bond vector before moving 0 and 1.
            drNext.Ev1Mv2(p1, p0);
            double px = p0.getX(0);
            double py = p0.getX(1);
            double fQ = 1 + b*drNext.squared();
            p0.PEa1Tv1(-omdth*fQ, drNext);
            if (a<0) {
                p0.setX(0, p0.getX(0) + srdt*py + sqdt*random.nextGaussian());
                p0.setX(1, p0.getX(1) + a*srdt*px + sqdt*random.nextGaussian());
            }
            else {
                p0.setX(0, p0.getX(0) + z*srdt*py + sqa*srdt*px + sqdt*random.nextGaussian());
                p0.setX(1, p0.getX(1) - sqa*srdt*py + sqdt*random.nextGaussian());
            }
            p0.setX(2, p0.getX(2) + sqdt*random.nextGaussian());
            center.PE(p0);

            for (int j=1; j<atoms.getAtomCount()-1; j++) {
                drPrev.E(drNext);
                p0 = p1;
                p1 = atoms.getAtom(j+1).getPosition();
                drNext.Ev1Mv2(p1, p0);
                px = p0.getX(0);
                py = p0.getX(1);

                fQ = 1 + b*drPrev.squared();
                p0.PEa1Tv1(fQ*omdth, drPrev);
                fQ = 1 + b*drNext.squared();
                p0.PEa1Tv1(-fQ*omdth, drNext);
                if (a<0) {
                    p0.setX(0, p0.getX(0) + srdt*py + sqdt*random.nextGaussian());
                    p0.setX(1, p0.getX(1) + a*srdt*px + sqdt*random.nextGaussian());
                }
                else {
                    p0.setX(0, p0.getX(0) + z*srdt*py + sqa*srdt*px + sqdt*random.nextGaussian());
                    p0.setX(1, p0.getX(1) - sqa*srdt*py + sqdt*random.nextGaussian());
                }
                p0.setX(2, p0.getX(2) + sqdt*random.nextGaussian());
                center.PE(p0);
            }

            drPrev.E(drNext);
            p0 = p1;
            px = p0.getX(0);
            py = p0.getX(1);

            fQ = 1 + b*drPrev.squared();
            p0.PEa1Tv1(fQ*omdth, drPrev);
            if (a<0) {
                p0.setX(0, p0.getX(0) + srdt*py + sqdt*random.nextGaussian());
                p0.setX(1, p0.getX(1) + a*srdt*px + sqdt*random.nextGaussian());
            }
            else {
                p0.setX(0, p0.getX(0) + z*srdt*py + sqa*srdt*px + sqdt*random.nextGaussian());
                p0.setX(1, p0.getX(1) - sqa*srdt*py + sqdt*random.nextGaussian());
            }
            p0.setX(2, p0.getX(2) + sqdt*random.nextGaussian());

            // maintain center at 0
            center.PE(p0);
            center.TE(-1.0/atoms.getAtomCount());

            for (int j=0; j<atoms.getAtomCount(); j++) {
                atoms.getAtom(j).getPosition().PE(center);
            }
        }
        
    }

    public void setTimeStep(double newTimeStep) {
        super.setTimeStep(newTimeStep);
        omdth = -0.25*timeStep;
        sqdt = Math.sqrt(0.5*timeStep);
    }

    public void setShearRateNumber(double newShearRate) {
        shearRate = newShearRate;
    }
    
    public double getShearRate() {
        if (a<0) return shearRate;
        double sqa = Math.sqrt(a);
        double sr = shearRate;
        if (sqa > 1.0/160.0) {
            sr /= 162*sqa;
        }
        return sr;
    }

    public double getShearRateNumber() {
        return shearRate;
    }

    public double getA() {
        return a;
    }

    public void setA(double newA) {
        a = newA;
    }

    public double getB() {
        return b;
    }

    public void setB(double newB) {
        b = newB;
    }

    private static final long serialVersionUID = 1L;
    protected final IVectorMutable drPrev, drNext, center;
    protected double omdth, sqdt;
    protected double shearRate;
    protected double a, b;
}
