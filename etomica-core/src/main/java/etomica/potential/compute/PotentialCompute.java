package etomica.potential.compute;

import etomica.atom.IAtom;
import etomica.molecule.IMolecule;
import etomica.space.Vector;

public interface PotentialCompute {
    void init();

    Vector[] getForces();

    double getLastVirial();

    double getOldEnergy();

    void updateAtom(IAtom atom);

    double computeAll(boolean doForces);

    double computeOneOld(IAtom iAtom);

    double computeOneOldMolecule(IMolecule molecule);

    double computeOne(IAtom iAtom);

    double computeOneMolecule(IMolecule molecule);

    void processAtomU(double fac);

    static PotentialCompute aggregate(PotentialCompute... computes) {
        return new PotentialComputeAggregate(computes);
    }
}
