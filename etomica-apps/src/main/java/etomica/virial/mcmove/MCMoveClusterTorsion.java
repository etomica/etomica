package etomica.virial.mcmove;

import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIterator;
import etomica.box.Box;
import etomica.data.histogram.HistogramSimple;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.integrator.mcmove.MCMoveStep;
import etomica.math.DoubleRange;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculePositionCOMPBC;
import etomica.potential.IPotentialBondTorsion;
import etomica.potential.PotentialMaster;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.util.collections.IntArrayList;
import etomica.util.random.IRandom;
import etomica.virial.BoxCluster;

import java.util.List;

public class MCMoveClusterTorsion extends MCMoveBoxStep {

    private final PotentialCompute potential;
    protected final IRandom random;
    protected final Space space;
    protected ISpecies species;
    protected double dt = 0;
    protected int a,b;
    protected Vector[][] position = null;
    protected final IntArrayList [] bonding;
    double uOld = 0;
    double uNew = 0;
    double wOld = 0;
    double wNew = 0;
    int [] indices = new int[3];
    int [] modified;
    int modifiedIndex = 0;
    HistogramSimple histogramSimple = new HistogramSimple(100, new DoubleRange(0, Math.PI));

    public MCMoveClusterTorsion(PotentialCompute potentialCompute, Space space, IntArrayList[] bonding, IRandom random, double stepSize) {
        super(null);
        this.potential = potentialCompute;
        this.space = space;
        this.random = random;
        this.bonding = bonding;
        this.stepSize = stepSize;
        setStepSizeMin(stepSize);
        setStepSizeMax(2*Math.PI);
        modified = new int[bonding.length];
    }

    public void setBox(Box p) {
        super.setBox(p);
        position = new Vector[p.getMoleculeList().size()][p.getMoleculeList().get(0).getChildList().size()];
    }

    @Override
    public boolean doTrial() {
        uOld = potential.computeAll(false);
        wOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        IMoleculeList moleculeList = box.getMoleculeList();
        for(int i = 0; i<moleculeList.size(); i++) {
            Vector com = space.makeVector();
            if (species != null && moleculeList.get(i).getType() != species) {
                continue;
            }
            IMolecule molecule = moleculeList.get(i);
            for(int j = 0; j < molecule.getChildList().size(); j++) {
                position[i][j] = space.makeVector();
                position[i][j].E(molecule.getChildList().getAtoms().get(j).getPosition());
            }
            com.E(MoleculePositionCOMPBC.com(box.getBoundary(), molecule));
            modifiedIndex = 0;
            dt = stepSize * (random.nextDouble() - 0.5);
            do{
                a = random.nextInt(bonding.length);
            }while (bonding[a].size() < 2);
            modified[modifiedIndex]=a;
            ++modifiedIndex;
            do{
                b = bonding[a].getInt(random.nextInt(bonding[a].size()));
            }while(bonding[b].size() < 2);
            modified[modifiedIndex] = b;
            ++modifiedIndex;
            List<IAtom> atoms = molecule.getChildList().getAtoms();
            for (int j = 0; j < bonding[a].size(); j++ ) {
                if (bonding[bonding[a].getInt(j)].size() > 1 && bonding[a].getInt(j) != b) {
                    for (int k = 0; k < bonding[b].size(); k++) {
                        if (bonding[bonding[b].getInt(k)].size() > 1 && bonding[b].getInt(k) != a) {
                            histogramSimple.addValue(computeTorsion(box.getBoundary(), atoms.get(bonding[a].getInt(j)), atoms.get(a), atoms.get(b), atoms.get(bonding[b].getInt(k))));
                        }
                    }
                }
            }
            Vector axis = space.makeVector();
            axis.Ev1Mv2(molecule.getChildList().getAtoms().get(a).getPosition(), molecule.getChildList().getAtoms().get(b).getPosition());
            axis.normalize();
            RotationTensor3D rotationTensor = new RotationTensor3D();
            rotationTensor.setRotationAxis(axis, dt/2);
            transformBondedAtoms(rotationTensor, position[i][a], a, molecule);
            rotationTensor.invert();
            transformBondedAtoms(rotationTensor, position[i][a], b, molecule);
            com.ME(MoleculePositionCOMPBC.com(box.getBoundary(), molecule));
            for(int j = 0; j < molecule.getChildList().size(); j++) {
                molecule.getChildList().getAtoms().get(j).getPosition().PE(com);
            }
        }
        ((BoxCluster)box).trialNotify();
        uNew = potential.computeAll(false);
        wNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        return true;
    }

    protected void transformBondedAtoms(RotationTensor3D rotationTensor3D, Vector pointOnAxis, int index, IMolecule molecule){
        for(int k = 0; k < bonding[index].size(); k++){
            boolean rotated = false;
            for (int l = 0; l < modifiedIndex; l++) {
                if (bonding[index].getInt(k) == modified[l]) {
                    rotated = true;
                    break;
                }
            }
            if (!rotated) {
                transform(rotationTensor3D, pointOnAxis, bonding[index].getInt(k), molecule);
                modified[modifiedIndex] = bonding[index].getInt(k);
                ++modifiedIndex;
                transformBondedAtoms(rotationTensor3D, pointOnAxis, bonding[index].getInt(k), molecule);
            }
        }
    }

    protected void transform(RotationTensor3D rotationTensor3D, Vector pointOnAxis, int index, IMolecule molecule) {
        //System.out.println("index " + index + " pointOnAxis " + pointOnAxis);
        Vector r = space.makeVector();
        r.Ev1Mv2(molecule.getChildList().getAtoms().get(index).getPosition(), pointOnAxis);
        rotationTensor3D.transform(r);
        r.PE(pointOnAxis);
        molecule.getChildList().getAtoms().get(index).getPosition().E(r);
    }

    private static double computeTorsion(Boundary boundary, IAtom iAtom, IAtom jAtom, IAtom kAtom, IAtom lAtom) {
        Vector ri = iAtom.getPosition();
        Vector rj = jAtom.getPosition();
        Vector rk = kAtom.getPosition();
        Vector rl = lAtom.getPosition();
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
        double theta = Math.acos(vji.dot(vkl) * vji_vkl);

        return theta;
    }

    public void printHistogram(){
        double [] histogram = histogramSimple.getHistogram();
        for (int i = 0; i < histogram.length; i++) {
            System.out.println(histogramSimple.xValues()[i] + "\t" + histogram[i]);
        }
    }

    @Override
    public double getChi(double temperature) {
        return (wOld == 0 ? 1 : wNew / wOld) * Math.exp(-(uNew - uOld) / temperature);
    }

    @Override
    public void acceptNotify() {
        ((BoxCluster)box).acceptNotify();
    }

    @Override
    public void rejectNotify() {
        IMoleculeList moleculeList = box.getMoleculeList();
        for(int i = 0; i<box.getMoleculeList().size(); i++) {
            for(int j = 0; j < moleculeList.get(i).getChildList().size(); j++) {
                moleculeList.get(i).getChildList().getAtoms().get(j).getPosition().E(position[i][j]);
            }
        }
        ((BoxCluster)box).rejectNotify();
    }

    @Override
    public AtomIterator affectedAtoms() {
        return null;
    }

    @Override
    public double energyChange() {
        return 0;
    }
}
