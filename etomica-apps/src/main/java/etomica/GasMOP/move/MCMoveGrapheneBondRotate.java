package etomica.GasMOP.move;

import etomica.atom.AtomSource;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.molecule.IMolecule;
import etomica.potential.compute.PotentialCompute;
import etomica.space.RotationTensor;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.collections.IntArrayList;
import etomica.util.random.IRandom;

import java.util.ArrayList;
import java.util.List;

public class MCMoveGrapheneBondRotate extends MCMoveBoxStep {
    protected final PotentialCompute potentialCompute;
    protected final List<Vector> oldPositions = new ArrayList<>();
    protected final IRandom random;
    protected IMolecule molecule;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected int maxIntArrSize;
    protected IntArrayList atomList;
    protected AtomSource atomSource;
    protected Space space;
    protected final Vector r0;
    protected final RotationTensor rotationTensor;
    protected final IntArrayList[] intArrayListMoeties;
    public MCMoveGrapheneBondRotate(IRandom random, PotentialCompute potentialCompute, Box box, IntArrayList[] intArrayListMoeties, IMolecule molecule) {
        super();
        this.potentialCompute = potentialCompute;
        this.random = random;
        this.space = box.getSpace();
        this.intArrayListMoeties = intArrayListMoeties;
        this.molecule = molecule;
        setStepSizeMax(Math.PI);
        setStepSizeMin(0.0);
        setStepSize(Math.PI / 4);
        perParticleFrequency = true;
        setBox(box);
        r0 = space.makeVector();
        rotationTensor = space.makeRotationTensor();
        maxIntArrSize = intArrayListMoeties.length;
    }

    @Override
    public boolean doTrial() {
        uOld = potentialCompute.computeAll(false);
        if (intArrayListMoeties.length==0) return false;
        atomList = intArrayListMoeties[random.nextInt(maxIntArrSize)];
        while (oldPositions.size() < atomList.size()) {
            oldPositions.add(space.makeVector());
        }
        for (int i = 0; i < atomList.size(); i++){
            Vector atomA = molecule.getChildList().get(atomList.getInt(i)).getPosition();
            oldPositions.get(i).E(atomA);
        }

        r0.E(molecule.getChildList().get(atomList.getInt(0)).getPosition());

        doTransform();
        return false;
    }

    @Override
    public double getChi(double temperature) {
        uNew = potentialCompute.computeOneMolecule(molecule);
//        System.out.println("rotate "+molecule.getIndex()+" "+uOld+" => "+uNew);
        return Math.exp(-(uNew - uOld) / temperature);
    }


    public double energyChange() {
        return uNew - uOld;
    }

    public void acceptNotify() {
//        System.out.println("accepted");
        potentialCompute.processAtomU(1);
        // put it back, then compute old contributions to energy
        IAtom atom;
        for (int i = 1; i < atomList.size(); i++){
            Vector atomA = molecule.getChildList().get(atomList.getInt(i)).getPosition();
            atom = molecule.getChildList().get(atomList.getInt(i));
            atom.getPosition().E(oldPositions.get(atom.getIndex()));
            potentialCompute.updateAtom(atom);
            oldPositions.get(i).E(atomA);
        }
        potentialCompute.computeOneMolecule(molecule);
        potentialCompute.processAtomU(-1);
        doTransform();
    }

    public void rejectNotify() {
        IAtom atom;
        for (int i = 1; i < atomList.size(); i++){
            Vector atomA = molecule.getChildList().get(atomList.getInt(i)).getPosition();
            atom = molecule.getChildList().get(atomList.getInt(i));
            atom.getPosition().E(oldPositions.get(atom.getIndex()));
            potentialCompute.updateAtom(atom);
            oldPositions.get(i).E(atomA);
        }
    }

    public void setBox(Box p) {
        super.setBox(p);
    }

    protected void doTransform() {
        IAtom a;
        for (int i = 1; i < atomList.size(); i++){
            a = molecule.getChildList().get(atomList.getInt(i));
            Vector r = molecule.getChildList().get(atomList.getInt(i)).getPosition();
            r.ME(r0);
            box.getBoundary().nearestImage(r);
            rotationTensor.transform(r);
            r.PE(r0);
            r.PE(box.getBoundary().centralImage(r));
            potentialCompute.updateAtom(a);
        }
    }
}
