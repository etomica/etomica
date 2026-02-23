package etomica.GasMOP.move;


import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculeSource;
import etomica.molecule.MoleculeSourceRandomMolecule;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.RotationTensor3D;
import etomica.util.collections.IntArrayList;
import etomica.util.random.IRandom;

import java.util.ArrayList;
import java.util.List;

public class MCMoveGrapheneTorsion extends MCMoveBoxStep {
    protected final PotentialCompute potentialCompute;
    protected final List<Vector> oldPositions = new ArrayList<>();
    protected final IRandom random;
    protected IMolecule molecule;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected MoleculeSource moleculeSource;
    protected Space space;
    protected IntArrayList[] torsions;
    protected int arrNum ;
    protected IntArrayList arrTorsion;
    protected double PI_2, phi;
    protected final Vector r0;
    protected Vector v1, v2, v3, v4;
    protected final RotationTensor3D rotationTensor;
    public MCMoveGrapheneTorsion(IRandom random, PotentialCompute potentialCompute, Box box, IntArrayList[] intArrayLists){
        super();
        this.potentialCompute = potentialCompute;
        this.random = random;
        this.space = box.getSpace();
        this.torsions = intArrayLists;
        MoleculeSourceRandomMolecule source = new MoleculeSourceRandomMolecule();
        source.setRandomNumberGenerator(random);
        moleculeSource = source;
        setBox(box);
        PI_2 = 2 * 3.14;
        r0 = space.makeVector();
        rotationTensor = (RotationTensor3D) space.makeRotationTensor();
    }

    public boolean doTrial(){
        molecule = moleculeSource.getMolecule();
        if (molecule == null) return false;
        uOld = potentialCompute.computeOneOldMolecule(molecule);
        if (uOld > 1e10) {
            throw new RuntimeException("molecule " + molecule + " in box " + box + " has an overlap ("+uOld+")");
        }
        while (oldPositions.size() < molecule.getChildList().size()) {
            oldPositions.add(space.makeVector());
        }
        for (IAtom a : molecule.getChildList()) {
            oldPositions.get(a.getIndex()).E(a.getPosition());
        }
        phi = PI_2 * random.nextDouble();
        arrNum =  random.nextInt(torsions.length);
        arrTorsion = torsions[arrNum];
        v2 = molecule.getChildList().get(arrTorsion.getInt(2)).getPosition();
        v3 = molecule.getChildList().get(arrTorsion.getInt(3)).getPosition();
        v2.ME(v3);
        box.getBoundary().nearestImage(v2);

        double axis2 = v2.squared();
        if (axis2 == 0.0) {
            return false;
        }
        v2.TE(1.0 / Math.sqrt(axis2));

        double dPhi = (2 * random.nextDouble() - 1.0) * stepSize;

        rotationTensor.setRotationAxis(v2, dPhi);
        doTransform();
        return true;
    }

    public void doTransform(){
        v2.ME(v3);
        box.getBoundary().nearestImage(v2);

        //? dont know
    }

    public double getChi(double temperature){
        uNew = potentialCompute.computeOneMolecule(molecule);
        return Math.exp(-(uNew - uOld) / temperature);
    }


    public double energyChange() {
        return uNew - uOld;
    }

    public void acceptNotify() {
//        System.out.println("accepted");
        potentialCompute.processAtomU(1);
        // put it back, then compute old contributions to energy
        molecule.getChildList().forEach(atom -> {
            atom.getPosition().E(oldPositions.get(atom.getIndex()));
            potentialCompute.updateAtom(atom);
        });
        potentialCompute.computeOneMolecule(molecule);
        potentialCompute.processAtomU(-1);
        doTransform();
    }

    public void rejectNotify() {
//        System.out.println("rejected");
      /*  molecule.getChildList().forEach(atom -> {
            atom.getPosition().E(oldPositions.get(atom.getIndex()));
            potentialCompute.updateAtom(atom);
        });*/
    }


    public MoleculeSource getMoleculeSource() {
        return moleculeSource;
    }

    public static void main(String[] args) {

    }
}

