package etomica.potential.compute;

import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.integrator.IntegratorListener;
import etomica.molecule.IMolecule;
import etomica.space.Vector;

public interface PotentialCompute {
    void init();

    default boolean needForcesForVirial() {
        return false;
    }

    ;

    Vector[] getForces();

    double getLastVirial();

    double getLastEnergy();

    void updateAtom(IAtom atom);

    default double computeAll(boolean doForces) {
        return computeAll(doForces, null);
    }

    double computeAll(boolean doForces, PotentialCallback pc);

    double computeOneOld(IAtom iAtom);

    double computeOneOldMolecule(IMolecule molecule);

    double computeOne(IAtom iAtom);

    default double computeOneMolecule(IMolecule molecule) {
        return computeManyAtoms(((AtomArrayList) molecule.getChildList()).toArray());
    }

    double computeManyAtoms(IAtom... atoms);

    void processAtomU(double fac);

    IntegratorListener makeIntegratorListener();

    static PotentialCompute aggregate(PotentialCompute... computes) {
        return new PotentialComputeAggregate(computes);
    }
}
