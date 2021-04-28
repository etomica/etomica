package etomica.virial;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.molecule.IMoleculeList;
import etomica.potential.Potential2Soft;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.SpeciesManager;

public class PotentialComputeIntramolecular implements PotentialCompute {
    public static boolean doDebug = false;
    protected final Space space;
    protected final Box box;
    protected final boolean isPureAtoms;
    protected final PotentialMasterBonding.FullBondingInfo bondingInfo;
    protected final Potential2Soft[][] atomPotentials;
    protected double u;
    protected Vector [] forces = null;

    public PotentialComputeIntramolecular(Space space, Box box, SpeciesManager sm, PotentialMasterBonding.FullBondingInfo bondingInfo, Potential2Soft[][] potential2Softs) {
        this.bondingInfo = bondingInfo;
        this.space = space;
        this.box = box;
        // the species we apply to might be purely atomic.  so long as our
        // BondingInfo is nonBonding, it won't matter
        isPureAtoms = sm.isPureAtoms();
        atomPotentials = potential2Softs;
    }
    @Override
    public void init() {

    }

    @Override
    public Vector[] getForces() {
        return forces;
    }

    @Override
    public double getLastVirial() {
        return 0;
    }

    @Override
    public double getLastEnergy() {
        return u;
    }

    @Override
    public void updateAtom(IAtom atom) {

    }

    @Override
    public double computeAll(boolean doForces, PotentialCallback pc) {
        IMoleculeList molecules = box.getMoleculeList();
        forces = space.makeVectorArray(box.getLeafList().size());

        u = 0;
        for (int k = 0; k < molecules.size(); k++){
            IAtomList atoms = molecules.get(k).getChildList();
            for (int i = 0; i < atoms.size(); i++) {
                IAtom a0 = atoms.get(i);
                Potential2Soft[] p0 = atomPotentials[a0.getType().getIndex()];
                for (int j = i + 1; j < atoms.size(); j++) {
                    IAtom a1 = atoms.get(j);
                    if (bondingInfo.skipBondedPair(isPureAtoms, a0, a1)) continue;
                    Potential2Soft p2 = p0[a1.getType().getIndex()];
                    if (p2 == null) continue;
                    Vector dr = space.makeVector();
                    dr.Ev1Mv2(a1.getPosition(), a0.getPosition());
                    box.getBoundary().nearestImage(dr);
                    double[] u012 = new double[3];
                    double r2 = dr.squared();
                    p2.u012add(r2, u012);
                    double uij = u012[0];
                    if (doDebug)
                        System.out.printf("i %d, j%d, uij %f \n", i, j, uij);
                    if (uij == 0) continue;
                    double duij = u012[1];
                    //System.out.printf("Force: %f", duij);
                    //System.out.println();
                    if (doForces) {
                        dr.TE(duij / r2);
                        forces[a0.getLeafIndex()].PE(dr);
                        forces[a1.getLeafIndex()].ME(dr);
                    }
                    u += uij;
                }
            }
        }
        //System.out.printf("PIntra %f \n", u);
        return u;
    }

    @Override
    public double computeOneOld(IAtom iAtom) {
        return 0;
    }

    @Override
    public double computeOne(IAtom iAtom) {
        return 0;
    }

    @Override
    public double computeManyAtomsOld(IAtom... atoms) {
        return 0;
    }

    @Override
    public double computeManyAtoms(IAtom... atoms) {
        return 0;
    }

    @Override
    public void processAtomU(double fac) {

    }

    @Override
    public IntegratorListener makeIntegratorListener() {
        return new IntegratorListener() {
            @Override
            public void integratorInitialized(IntegratorEvent e) {

            }
        };
    }
}
