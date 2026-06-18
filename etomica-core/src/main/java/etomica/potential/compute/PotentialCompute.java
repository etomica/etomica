package etomica.potential.compute;

import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorListener;
import etomica.molecule.CenterOfMass;
import etomica.molecule.IMolecule;
import etomica.space.Vector;

public interface PotentialCompute {
    void init();

    default boolean needForcesForVirial() {
        return false;
    }

    default Vector[] getTorques() {return null;}

    Vector[] getForces();

    double getLastVirial();

    static double computeVirialIntramolecular(Vector[] forces, Box box) {
        double virialIntra = 0;
        for (IMolecule molecule : box.getMoleculeList()) {
            IAtomList atoms = molecule.getChildList();
            if (atoms.size() == 1) continue;
            Vector com = CenterOfMass.position(box, molecule);
            for (IAtom atom : atoms) {
                Vector ri = atom.getPosition();
                Vector dr = box.getSpace().makeVector();
                dr.Ev1Mv2(ri, com);
                box.getBoundary().nearestImage(dr);
                virialIntra += forces[atom.getLeafIndex()].dot(dr);
            }
        }
        return virialIntra;
    }

    double getLastEnergy();

    void updateAtom(IAtom atom);

    default double computeAll(boolean doForces) {
        return computeAll(doForces, null);
    }

    double computeAll(boolean doForces, PotentialCallback pc);

    double computeOneOld(IAtom iAtom);

    default double computeOneOldMolecule(IMolecule molecule) {
        return computeManyAtomsOld(((AtomArrayList) molecule.getChildList()).toArray());
    }

    double computeOne(IAtom iAtom);

    default double computeOneMolecule(IMolecule molecule) {
        return computeManyAtoms(((AtomArrayList) molecule.getChildList()).toArray());
    }

    double computeManyAtomsOld(IAtom... atoms);

    double computeManyAtoms(IAtom... atoms);

    void processAtomU(double fac);

    IntegratorListener makeIntegratorListener();

    static PotentialCompute aggregate(PotentialCompute... computes) {
        return new PotentialComputeAggregate(computes);
    }
}
