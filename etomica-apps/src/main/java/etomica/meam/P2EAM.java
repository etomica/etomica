package etomica.meam;

import etomica.api.*;
import etomica.atom.iterator.IteratorDirective;
import etomica.potential.Potential2;
import etomica.potential.PotentialCalculationEnergySum;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialSoft;
import etomica.space.ISpace;
import etomica.space.Tensor;

/**
 * EAM (Embedded Atom Method) potential that pretends to be a pair potential.
 * Since it is only pretending, the methods must be called in a somewhat
 * different manner.  The different structure should result in higher
 * performance, as compared to {@link PotentialEAM}.
 * <p>
 * To compute the energy, the P1EAM inner class should also be added to the
 * potential master.  The energy method should be called for all interacting
 * pairs and then summed with the result of calling energy1.
 * <p>
 * To compute the gradient, the energy method should be called for all pairs.
 * Then call prepForGradient.  Then the gradient method should be called for
 * all pairs and summed appropriately.  To speed things up, disableEnergy may
 * be called.
 *
 * @author Andrew Schultz
 */
public class P2EAM extends Potential2 implements PotentialSoft {
    
    protected double n2, m2, eps, a, a2, Ceps, rc12, rc22;
    protected IBoundary boundary;
    protected final IVectorMutable dr;
    protected IVectorMutable[] gradient;
    protected IVectorMutable rhograd;
    protected double[] rho;
    protected boolean energyDisabled = false;

    public P2EAM(ISpace space, double n, double m, double eps, double a, double C, double rc1, double rc2) {
        super(space);

        n2 = 0.5 * n;
        m2 = 0.5 * m;
        this.eps = eps;
        this.a = a;
        a2 = a * a;
        Ceps = C*eps;
        rc12 = rc1 * rc1;
        rc22 = rc2 * rc2;
        dr = space.makeVector();
        gradient = new IVectorMutable[2];
        gradient[0] = space.makeVector();
        gradient[1] = space.makeVector();
        rhograd = space.makeVector();
        rho = new double[0];
    }

    public double getRange() {
        return Math.sqrt(rc12 > rc22 ? rc12 : rc22);
    }

    public void setCutoff(double rc1, double rc2) {
        rc12 = rc1 * rc1;
        rc22 = rc2 * rc2;
    }

    public void disableEnergy() {
        energyDisabled = true;
    }

    public void enableEnergy() {
        energyDisabled = false;
    }

    public void reset() {
        for (int i = 0; i < rho.length; i++) {
            rho[i] = 0;
        }
    }

    public void prepForGradient() {
        for (int i = 0; i < rho.length; i++) {
            rho[i] = 1 / Math.sqrt(rho[i]);
        }
    }

    public double energy(IAtomList atoms) {
        IVector pos0 = atoms.getAtom(0).getPosition();
        IVector pos1 = atoms.getAtom(1).getPosition();
        dr.Ev1Mv2(pos1, pos0);
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        double u = 0;
        if (!energyDisabled && r2 < rc12) {
            u = eps * Math.pow(a2 / r2, n2);
        }
        if (r2 < rc22) {
            double rhoi = Math.pow(a2 / r2, m2);
            rho[atoms.getAtom(0).getLeafIndex()] += rhoi;
            rho[atoms.getAtom(1).getLeafIndex()] += rhoi;
        }

        return u;
    }
    
    public double energy1() {
        double sum1 = 0;
        for (int i=0; i<rho.length; i++) {
            sum1 += -Ceps * Math.sqrt(rho[i]);
        }
        return sum1;
    }

    public void setBox(IBox box) {
        boundary = box.getBoundary();
        if (rho.length != box.getLeafList().getAtomCount()) {
            rho = new double[box.getLeafList().getAtomCount()];
        }
    }

    public double virial(IAtomList atoms) {
        throw new RuntimeException("implement me");
    }

    public IVector[] gradient(IAtomList atoms) {

        IVector pos0 = atoms.getAtom(0).getPosition();
        IVector pos1 = atoms.getAtom(1).getPosition();
        dr.Ev1Mv2(pos1, pos0);
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        if (r2 > rc12 && r2 > rc22) {
            gradient[0].E(0);
            gradient[1].E(0);
            return gradient;
        }
        double a2r2 = a2/r2;
        if (r2 < rc12) {
            double u = eps * Math.pow(a2r2, n2);
            double dvdr = -2 * n2 * u;
            gradient[1].Ea1Tv1(-dvdr, dr);
        }
        else {
            gradient[1].E(0);
        }

        if (r2 <= rc22) {
            double drhodr = -2 * m2 * Math.pow(a2r2, m2);
            rhograd.Ea1Tv1(drhodr * Ceps / 2, dr);
            double sqrtRhoDiff = rho[atoms.getAtom(1).getLeafIndex()] - rho[atoms.getAtom(0).getLeafIndex()];
            gradient[1].PEa1Tv1(sqrtRhoDiff, rhograd);
        }
        gradient[0].Ea1Tv1(-1, gradient[1]);
        return gradient;
    }

    public IVector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        throw new RuntimeException("not implemented.  use gradient");
    }
    
    /**
     * Creates and returns an integrator listener that makes the integrator's
     * force sum work properly by running an energy calculation first.
     *
     * @param potentialMaster
     * @param box
     * @return the listener
     */
    public IIntegratorListenerMD makeIntegratorListener(PotentialMaster potentialMaster, IBox box) {
        final P2EAM p2 = this;
        return new IIntegratorListenerMD() {
            PotentialCalculationEnergySum pcEnergy = new PotentialCalculationEnergySum();
            IteratorDirective id = new IteratorDirective();
            @Override
            public void integratorForcePrecomputed(IIntegratorEvent e) {
                p2.disableEnergy();
                p2.reset();
                potentialMaster.calculate(box, id, pcEnergy);
                p2.enableEnergy();
            }
        
            @Override
            public void integratorForceComputed(IIntegratorEvent e) {}
        
            @Override
            public void integratorInitialized(IIntegratorEvent e) {}
        
            @Override
            public void integratorStepStarted(IIntegratorEvent e) {}
        
            @Override
            public void integratorStepFinished(IIntegratorEvent e) {}
        };
    }
}
