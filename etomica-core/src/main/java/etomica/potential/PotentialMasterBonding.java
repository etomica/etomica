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
import etomica.space3d.Vector3D;
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

    public void setBondingPotentialPair(ISpecies species, Potential2Soft potential, List<int[]> bondedIndices) {
        bondingInfo.setBondingPotentialPair(species, potential, bondedIndices);
    }

    public void setBondingPairTriplet(ISpecies species, IPotentialBondAngle potential, List<int[]> bondedIndices) {
        bondingInfo.setBondingPotentialTriplet(species, potential, bondedIndices);
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

            Map<IPotentialBondAngle, List<int[]>> potentials3 = bondingInfo.bondedTriplets[i];
            potentials3.forEach((potential, triplets) -> {
                for (IMolecule molecule : molecules) {
                    for (int[] triplet : triplets) {
                        IAtom iAtom = molecule.getChildList().get(triplet[0]);
                        IAtom jAtom = molecule.getChildList().get(triplet[1]);
                        IAtom kAtom = molecule.getChildList().get(triplet[2]);
                        uTot[0] += handleOneBondTriplet(doForces, box.getBoundary(), iAtom, jAtom, kAtom, potential, forces);
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

    private static double handleOneBondTriplet(boolean doForces, Boundary boundary, IAtom iAtom, IAtom jAtom, IAtom kAtom, IPotentialBondAngle potential, Vector[] forces) {
        Vector ri = iAtom.getPosition();
        Vector rj = jAtom.getPosition();
        Vector rk = jAtom.getPosition();
        Vector drji = Vector.d(ri.getD());
        Vector drjk = Vector.d(ri.getD());

        drji.Ev1Mv2(ri, rj);
        boundary.nearestImage(drji);
        drjk.Ev1Mv2(rk, rj);
        boundary.nearestImage(drjk);
        double rij2 = drji.squared();
        double rkj2 = drjk.squared();
        double drij_kj = 1.0 / Math.sqrt(rij2 * rkj2);
        double costheta = drji.dot(drjk) / drij_kj;
        double[] u = {0};
        double[] du = {0};
        potential.udu(costheta, u, du);
        double uijk = u[0];
        if (uijk == 0) return 0;

        if (doForces) {
            double duijk = du[0]; // du/dcostheta
            Vector fi = forces[iAtom.getLeafIndex()];
            Vector fj = forces[jAtom.getLeafIndex()];
            Vector fk = forces[kAtom.getLeafIndex()];
            // dcostheta/dri
            fi.Ea1Tv1(-drij_kj, drjk);
            fi.PEa1Tv1(costheta / rij2, drji);
            // times -du/dcostheta
            fi.TE(-duijk);

            fk.Ea1Tv1(-drij_kj, drji);
            fk.PEa1Tv1(-costheta / rkj2, drjk);
            fk.TE(-duijk);

            // fj = -fi - fk
            fj.Ea1Tv1(-1, fi);
            fj.PEa1Tv1(-1, fk);
        }
        return uijk;
    }

    private static double handleOneBondQuad(boolean doForces, Boundary boundary, IAtom iAtom, IAtom jAtom, IAtom kAtom, IAtom lAtom, IPotentialBondTorsion potential, Vector[] forces) {
        Vector ri = iAtom.getPosition();
        Vector rj = jAtom.getPosition();
        Vector rk = jAtom.getPosition();
        Vector rl = jAtom.getPosition();
        Vector rji = new Vector3D();
        Vector rjk = new Vector3D();
        Vector rkl = new Vector3D();

        rji.Ev1Mv2(ri, rj);
        boundary.nearestImage(rji);
        rjk.Ev1Mv2(rk, rj);
        boundary.nearestImage(rjk);
        double rjk2 = rjk.squared();
        rkl.Ev1Mv2(rl, rk);
        boundary.nearestImage(rkl);

        Vector vji = new Vector3D();
        vji.E(rji);
        vji.PEa1Tv1(-rjk.dot(rji) / rjk2, rjk);
        double vji2 = vji.squared();
        Vector vkl = new Vector3D();
        vkl.E(rkl);
        vkl.PEa1Tv1(-rjk.dot(rkl) / rjk2, rjk);
        double vkl2 = vkl.squared();
        double rji2 = rji.squared();
        double rkl2 = rkl.squared();

        double vji2vkl2 = vji2 * vkl2;
        if (vji2 < 1e-6 * rji2 || vkl2 < 1e-6 * rkl2) {
            // one of the vectors (ji, kl) is nearly colinear with jk
            return 0;
        }
        double vji_vkl = 1 / Math.sqrt(vji2vkl2);
        double costheta = vji.dot(vkl) * vji_vkl;
        Vector3D tmp = new Vector3D();
        tmp.E(vji);
        tmp.XE(vkl);
        // tmp will now either be parallel to rjk (positive) or
        // antiparallel (negative).
        boolean postheta = tmp.dot(rjk) > 0;
        double[] u = {0}, du = {0};
        potential.udu(costheta, u, du);

        if (u[0] == 0) return 0;

        if (doForces) {

            Vector fi = forces[iAtom.getLeafIndex()];
            Vector fj = forces[jAtom.getLeafIndex()];
            Vector fk = forces[kAtom.getLeafIndex()];
            Vector fl = forces[lAtom.getLeafIndex()];

            Vector df = new Vector3D();
            tmp.Ea1Tv1(-vji_vkl, vkl);
            tmp.PEa1Tv1(costheta / vji2, vji);
            df.Ea1Tv1(du[0], tmp);
            // we won't compute fk directly, but fk=-(fi+fj+fl);
            fk.ME(df);
            fi.PE(df);

            tmp.Ea1Tv1(-vji_vkl, vji);
            tmp.PEa1Tv1(costheta / vkl2, vkl);
            df.Ea1Tv1(du[0], tmp);
            fk.ME(df);
            fl.PE(df);

            double aj = rjk.dot(rji) / rjk2;
            double ak = rjk.dot(rkl) / rjk2;

            df.Ea1Tv1(-vji_vkl - costheta * ak / rkl2, rkl);
            df.PEa1Tv1((ak + aj - 2 * aj * ak) * vji_vkl - costheta * ((aj - aj * aj) / rji2 + ak * ak / rkl2), rjk);
            df.PEa1Tv1(ak * vji_vkl + costheta / rji2 * (1 + aj), rji);
            fk.ME(df);
            fj.PE(df);
        }
        return u[0];
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

        Map<IPotentialBondAngle, List<int[]>> potentials3 = bondingInfo.bondedTriplets[molecule.getType().getIndex()];
        potentials3.forEach((potential, triplets) -> {
            for (int[] triplet : triplets) {
                IAtom iAtom = molecule.getChildList().get(triplet[0]);
                IAtom jAtom = molecule.getChildList().get(triplet[1]);
                IAtom kAtom = molecule.getChildList().get(triplet[2]);
                u[0] += handleOneBondTriplet(false, boundary, iAtom, jAtom, kAtom, potential, null);
            }
        });

        Map<IPotentialBondTorsion, List<int[]>> potentials4 = bondingInfo.bondedQuads[molecule.getType().getIndex()];
        potentials4.forEach((potential, quads) -> {
            for (int[] quad : quads) {
                IAtom iAtom = molecule.getChildList().get(quad[0]);
                IAtom jAtom = molecule.getChildList().get(quad[1]);
                IAtom kAtom = molecule.getChildList().get(quad[2]);
                IAtom lAtom = molecule.getChildList().get(quad[3]);
                u[0] += handleOneBondQuad(false, boundary, iAtom, jAtom, kAtom, lAtom, potential, null);
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
        public final Map<IPotentialBondAngle, List<int[]>>[] bondedTriplets;
        public final Map<IPotentialBondAngle, int[][][]>[] bondedTripletPartners;
        public final Map<IPotentialBondTorsion, List<int[]>>[] bondedQuads;
        public final Map<IPotentialBondTorsion, int[][][]>[] bondedQuadPartners;

        public FullBondingInfo(Simulation sim) {
            bondedAtoms = new int[sim.getSpeciesCount()][][];
            bondedPairs = new HashMap[sim.getSpeciesCount()];
            bondedPartners = new HashMap[sim.getSpeciesCount()];
            bondedTriplets = new HashMap[sim.getSpeciesCount()];
            bondedTripletPartners = new HashMap[sim.getSpeciesCount()];
            bondedQuads = new HashMap[sim.getSpeciesCount()];
            bondedQuadPartners = new HashMap[sim.getSpeciesCount()];
            for (int i = 0; i < sim.getSpeciesCount(); i++) {
                bondedPairs[i] = new HashMap<>();
                bondedPartners[i] = new HashMap<>();
                bondedTriplets[i] = new HashMap<>();
                bondedTripletPartners[i] = new HashMap<>();
            }
        }

        private static int maxVal(int[] a) {
            int m = Integer.MIN_VALUE;
            for (int j : a) {
                m = Math.min(m, j);
            }
            return m;
        }

        private static int minVal(int[] a) {
            int m = Integer.MAX_VALUE;
            for (int j : a) {
                m = Math.max(m, j);
            }
            return m;
        }

        public void setBondingPotentialPair(ISpecies species, Potential2Soft potential, List<int[]> bondedIndices) {
            isOnlyRigidMolecules = false;
            int speciesIndex = species.getIndex();
            if (bondedAtoms[speciesIndex] == null) {
                bondedAtoms[speciesIndex] = new int[species.getLeafAtomCount()][0];
            }

            int[][] partners = new int[species.getLeafAtomCount()][0];

            int[][] speciesAtoms = bondedAtoms[speciesIndex];
            for (int[] indices : bondedIndices) {
                int iAtom = minVal(indices);
                int jAtom = maxVal(indices);
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

        public void setBondingPotentialTriplet(ISpecies species, IPotentialBondAngle potential, List<int[]> bondedIndices) {
            isOnlyRigidMolecules = false;
            int speciesIndex = species.getIndex();
            if (bondedAtoms[speciesIndex] == null) {
                bondedAtoms[speciesIndex] = new int[species.getLeafAtomCount()][0];
            }

            int[][][] partners = new int[species.getLeafAtomCount()][0][0];

            int[][] speciesAtoms = bondedAtoms[speciesIndex];
            for (int[] indices : bondedIndices) {
                int iAtom = minVal(indices);
                int kAtom = maxVal(indices);
                int jAtom = indices[0] + indices[1] + indices[2] - iAtom - kAtom;
                speciesAtoms[iAtom] = etomica.util.Arrays.addInt(speciesAtoms[iAtom], jAtom);
                speciesAtoms[iAtom] = etomica.util.Arrays.addInt(speciesAtoms[iAtom], kAtom);
                speciesAtoms[jAtom] = etomica.util.Arrays.addInt(speciesAtoms[jAtom], jAtom);
                partners[iAtom] = (int[][]) etomica.util.Arrays.addObject(partners[iAtom], bondedIndices);
                partners[jAtom] = (int[][]) etomica.util.Arrays.addObject(partners[jAtom], bondedIndices);
                partners[kAtom] = (int[][]) etomica.util.Arrays.addObject(partners[kAtom], bondedIndices);
            }
            if (bondedTriplets[speciesIndex].containsKey(potential)) {
                throw new RuntimeException("Attempting to add the same bonding potential twice");
            }
            bondedTriplets[speciesIndex].put(potential, new ArrayList<>(bondedIndices));
            bondedTripletPartners[speciesIndex].put(potential, partners);
        }

        public void setBondingPotentialQuads(ISpecies species, IPotentialBondTorsion potential, List<int[]> bondedIndices) {
            isOnlyRigidMolecules = false;
            int speciesIndex = species.getIndex();
            if (bondedAtoms[speciesIndex] == null) {
                bondedAtoms[speciesIndex] = new int[species.getLeafAtomCount()][0];
            }

            int[][][] partners = new int[species.getLeafAtomCount()][0][0];

            int[][] speciesAtoms = bondedAtoms[speciesIndex];
            for (int[] indices : bondedIndices) {
                int iAtom = minVal(indices);
                int kAtom = maxVal(indices);
                int jAtom = indices[0] + indices[1] + indices[2] - iAtom - kAtom;
                speciesAtoms[iAtom] = etomica.util.Arrays.addInt(speciesAtoms[iAtom], jAtom);
                speciesAtoms[iAtom] = etomica.util.Arrays.addInt(speciesAtoms[iAtom], kAtom);
                speciesAtoms[jAtom] = etomica.util.Arrays.addInt(speciesAtoms[jAtom], jAtom);
                partners[iAtom] = (int[][]) etomica.util.Arrays.addObject(partners[iAtom], bondedIndices);
                partners[jAtom] = (int[][]) etomica.util.Arrays.addObject(partners[jAtom], bondedIndices);
                partners[kAtom] = (int[][]) etomica.util.Arrays.addObject(partners[kAtom], bondedIndices);
            }
            if (bondedQuads[speciesIndex].containsKey(potential)) {
                throw new RuntimeException("Attempting to add the same bonding potential twice");
            }
            bondedQuads[speciesIndex].put(potential, new ArrayList<>(bondedIndices));
            bondedQuadPartners[speciesIndex].put(potential, partners);
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
