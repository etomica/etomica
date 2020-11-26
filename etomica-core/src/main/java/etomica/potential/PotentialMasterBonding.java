package etomica.potential;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.compute.PotentialCompute;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;

import java.util.*;

public class PotentialMasterBonding implements PotentialCompute {

    private final List<ISpecies> speciesList;
    protected Vector[] forces = new Vector[0];
    private final Box box;
    private final Space space;
    private final FullBondingInfo bondingInfo;

    public PotentialMasterBonding(Simulation sim, Box box) {
        this(sim, box, new FullBondingInfo(sim));
    }

    public PotentialMasterBonding(Simulation sim, Box box, FullBondingInfo bondingInfo) {
        speciesList = sim.getSpeciesList();
        this.box = box;
        this.space = box.getSpace();
        this.bondingInfo = bondingInfo;
    }

    public FullBondingInfo getBondingInfo() {
        return bondingInfo;
    }

    public void setBondingPotential(ISpecies species, Potential2Soft potential, List<int[]> bondedIndices) {
        bondingInfo.setBondingPotential(species, potential, bondedIndices);
    }

    private void zeroArrays(boolean doForces) {
        if (doForces) {
            int numAtoms = box.getLeafList().size();
            if (numAtoms > forces.length) {
                int oldLength = forces.length;
                forces = Arrays.copyOf(forces, numAtoms);
                for (int i = oldLength; i < numAtoms; i++) forces[i] = box.getSpace().makeVector();
            }
            for (int i = 0; i < numAtoms; i++) {
                forces[i].E(0);
            }
        }
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
    public double getOldEnergy() {
        return computeAll(false);
    }

    @Override
    public void updateAtom(IAtom atom) {

    }

    public double computeAll(boolean doForces) {
        zeroArrays(doForces);

        double[] uTot = {0};
        for (int i = 0; i < speciesList.size(); i++) {
            Map<Potential2Soft, List<int[]>> potentials = bondingInfo.bondedPairs[i];
            IMoleculeList molecules = box.getMoleculeList(speciesList.get(i));
            potentials.forEach((potential, pairs) -> {
                for (IMolecule molecule : molecules) {
                    for (int[] pair : pairs) {
                        IAtom iAtom = molecule.getChildList().get(pair[0]);
                        IAtom jAtom = molecule.getChildList().get(pair[1]);
                        uTot[0] += handleOneBondPair(doForces, box.getBoundary(), iAtom, jAtom, potential, forces);
                    }
                }
            });
        }
        return uTot[0];
    }

    public double computeOne(IAtom atom) {
        double[] uTot = {0};
        IMolecule parentMolecule = atom.getParentGroup();
        bondingInfo.bondedPartners[atom.getParentGroup().getType().getIndex()].forEach((potential, partners) -> {
            for (int partnerIdx : partners[atom.getIndex()]) {
                IAtom jAtom = parentMolecule.getChildList().get(partnerIdx);
                uTot[0] += handleComputeOneBonded(potential, atom, jAtom);
            }
        });
        return uTot[0];
    }

    private double handleComputeOneBonded(Potential2Soft pij, IAtom iAtom, IAtom jAtom) {
        Vector ri = iAtom.getPosition();
        Vector rj = jAtom.getPosition();
        Vector dr = space.makeVector();
        dr.Ev1Mv2(rj, ri);
        box.getBoundary().nearestImage(dr);
        return pij.u(dr.squared());
    }

    private static double handleOneBondPair(boolean doForces, Boundary boundary, IAtom iAtom, IAtom jAtom, Potential2Soft potential, Vector[] forces) {
        Vector ri = iAtom.getPosition();
        Vector rj = jAtom.getPosition();
        Vector dr = Vector.d(ri.getD());

        dr.Ev1Mv2(rj, ri);
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        double[] u = {0};
        double[] du = {0};
        potential.udu(r2, u, du);
        double uij = u[0];
        if (uij == 0) return 0;

        if (doForces) {
            double duij = du[0];
            dr.TE(duij / r2);
            forces[iAtom.getLeafIndex()].PE(dr);
            forces[jAtom.getLeafIndex()].ME(dr);
        }
        return uij;
    }

    public double computeOneMolecule(IMolecule molecule) {
        return computeOneMolecule(box.getBoundary(), molecule, bondingInfo);
    }

    public static double computeOneMolecule(Boundary boundary, IMolecule molecule, FullBondingInfo bondingInfo) {
        final double[] u = {0};
        Map<Potential2Soft, List<int[]>> potentials = bondingInfo.bondedPairs[molecule.getType().getIndex()];
        potentials.forEach((potential, pairs) -> {
            for (int[] pair : pairs) {
                IAtom iAtom = molecule.getChildList().get(pair[0]);
                IAtom jAtom = molecule.getChildList().get(pair[1]);
                u[0] += handleOneBondPair(false, boundary, iAtom, jAtom, potential, null);
            }
        });
        return u[0];
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

            @Override
            public void integratorStepStarted(IntegratorEvent e) {

            }

            @Override
            public void integratorStepFinished(IntegratorEvent e) {

            }
        };
    }

    public double computeOneOld(IAtom iAtom) {
        return computeOne(iAtom);
    }

    public double computeOneOldMolecule(IMolecule molecule) {
        if (!bondingInfo.isOnlyRigidMolecules) {
            throw new RuntimeException();
        }
        return computeOneMolecule(molecule);
    }

    public static class FullBondingInfo implements BondingInfo {
        public boolean isOnlyRigidMolecules = true;
        public final int[][][] bondedAtoms;
        public final Map<Potential2Soft, List<int[]>>[] bondedPairs;
        public final Map<Potential2Soft, int[][]>[] bondedPartners;

        public FullBondingInfo(Simulation sim) {
            bondedAtoms = new int[sim.getSpeciesCount()][][];
            bondedPairs = new HashMap[sim.getSpeciesCount()];
            bondedPartners = new HashMap[sim.getSpeciesCount()];
            for (int i = 0; i < bondedPairs.length; i++) {
                bondedPairs[i] = new HashMap<>();
                bondedPartners[i] = new HashMap<>();
            }
        }

        public void setBondingPotential(ISpecies species, Potential2Soft potential, List<int[]> bondedIndices) {
            isOnlyRigidMolecules = false;
            int speciesIndex = species.getIndex();
            if (bondedAtoms[speciesIndex] == null) {
                bondedAtoms[speciesIndex] = new int[species.getLeafAtomCount()][0];
            }

            int[][] partners = new int[species.getLeafAtomCount()][0];

            int[][] speciesAtoms = bondedAtoms[speciesIndex];
            for (int[] indices : bondedIndices) {
                int iAtom = Math.min(indices[0], indices[1]);
                int jAtom = Math.max(indices[0], indices[1]);
                speciesAtoms[iAtom] = etomica.util.Arrays.addInt(speciesAtoms[iAtom], jAtom);
                partners[iAtom] = etomica.util.Arrays.addInt(partners[iAtom], jAtom);
                partners[jAtom] = etomica.util.Arrays.addInt(partners[jAtom], iAtom);
            }
            if (bondedPairs[speciesIndex].containsKey(potential)) {
                throw new RuntimeException("Attempting to add the same bonding potential twice");
            }
            bondedPairs[speciesIndex].put(potential, new ArrayList<>(bondedIndices));
            bondedPartners[speciesIndex].put(potential, partners);
        }

        @Override
        public boolean skipBondedPair(boolean isPureAtoms, IAtom iAtom, IAtom jAtom) {
            if (!isPureAtoms && iAtom.getParentGroup() == jAtom.getParentGroup()) {
                // ensure i < j
                if (isOnlyRigidMolecules) {
                    return true;
                }
                if (iAtom.getLeafIndex() > jAtom.getLeafIndex()) {
                    IAtom tmp = iAtom;
                    iAtom = jAtom;
                    jAtom = tmp;
                }
                int species = iAtom.getParentGroup().getType().getIndex();
                int iChildIndex = iAtom.getIndex();
                int jChildIndex = jAtom.getIndex();
                int[] iBondedAtoms = bondedAtoms[species][iChildIndex];
                for (int iBondedAtom : iBondedAtoms) {
                    if (iBondedAtom == jChildIndex) return true;
                }
            }
            return false;
        }

        @Override
        public boolean isOnlyRigidMolecules() {
            return isOnlyRigidMolecules;
        }
    }
}
