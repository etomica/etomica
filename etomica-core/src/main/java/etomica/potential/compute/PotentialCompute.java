package etomica.potential.compute;

import etomica.atom.IAtom;
import etomica.integrator.IntegratorListener;
import etomica.molecule.IMolecule;
import etomica.space.Vector;

public interface PotentialCompute {
    void init();

    Vector[] getForces();

    double getLastVirial();

    double getOldEnergy();

    void updateAtom(IAtom atom);

    default double computeAll(boolean doForces) {
        return computeAll(doForces, null);
    }

    double computeAll(boolean doForces, PotentialCallback pc);

    double computeOneOld(IAtom iAtom);

    double computeOneOldMolecule(IMolecule molecule);

    double computeOne(IAtom iAtom);

    double computeOneMolecule(IMolecule molecule);

    void processAtomU(double fac);

    IntegratorListener makeIntegratorListener();

    static PotentialCompute aggregate(PotentialCompute... computes) {
        return new PotentialComputeAggregate(computes);
    }
}
