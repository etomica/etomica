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
import etomica.space3d.Vector3D;
import etomica.units.Degree;
import etomica.util.collections.IntArrayList;
import etomica.util.random.IRandom;

import java.util.ArrayList;
import java.util.List;

public class MCMoveGrapheneTorsion extends MCMoveBoxStep  {
    protected final PotentialCompute potentialCompute;
    protected final List<Vector> oldPositions = new ArrayList<>();
    public Box box;
    protected double stepSize = 1.0;
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
    protected Vector v0, v1, v2, v3, v01, v12;
    protected double r01_2, r12_2, dr01_12, theta;
    protected final RotationTensor3D rotationTensor;
    // keep angle same but change torsion. Explore torsions in GO
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
        v0 = molecule.getChildList().get(arrTorsion.getInt(0)).getPosition();
        v1 = molecule.getChildList().get(arrTorsion.getInt(1)).getPosition();
        v2 = molecule.getChildList().get(arrTorsion.getInt(2)).getPosition();
        v3 = molecule.getChildList().get(arrTorsion.getInt(3)).getPosition();

        v01.Ev1Mv2(v0, v1);
        v12.Ev1Mv2(v1, v2);
        r01_2 = v01.squared();
        r12_2 = v12.squared();
        dr01_12 = 1.0 / Math.sqrt(r12_2 * r01_2);
        theta = Math.acos( v01.dot(v12) * dr01_12);
        /*
        v2.ME(v3);
        box.getBoundary().nearestImage(v2);
        double axis2 = v2.squared();
        if (axis2 == 0.0) {
            return false;
        }
        v2.TE(1.0 / Math.sqrt(axis2));

        double dPhi = (2 * random.nextDouble() - 1.0) * stepSize;
        rotationTensor.setRotationAxis(v2, dPhi);*/
        doTransform();
        return true;
    }

    public void doTransform(){
        v2.ME(v3);
        box.getBoundary().nearestImage(v2);

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
        molecule.getChildList().forEach(atom -> {
            atom.getPosition().E(oldPositions.get(atom.getIndex()));
            potentialCompute.updateAtom(atom);
        });
    }


    public MoleculeSource getMoleculeSource() {
        return moleculeSource;
    }

    public static void main(String[] args) {
        Vector vec1 = new Vector3D(-7,5,3);
        Vector vec2 = new Vector3D(2,0,0);
        Vector vec3 = new Vector3D(1,5,0);
        Vector vec4 = new Vector3D(2,7,1);
        Vector v12 = new Vector3D();
        Vector v23 = new Vector3D();
        v12.E(vec1);
        v12.ME(vec2);
        System.out.println(v12);
        v23.E(vec2);
        v23.ME(vec3);
        System.out.println(v23);
        Vector v12Norm = new Vector3D();
        Vector v23Norm = new Vector3D();
        v12Norm.E(v12);
        v23Norm.E(v23);
        v12Norm.normalize();
        v23Norm.normalize();
        double cosTheta = v12Norm.dot(v23Norm);
        double theta = Math.acos(cosTheta);
        System.out.println(theta);

        v23 = new Vector3D();
        Vector v34 = new Vector3D();
        v23.E(vec2);
        v23.ME(vec3);
        v34.E(vec3);
        v34.ME(vec4);

        Vector n1Cross = new Vector3D();
        Vector n2Cross = new Vector3D();
        n1Cross.E(v12);
        n1Cross.XE(v23);

        n2Cross.E(v23);
        n2Cross.XE(v34);

        Vector v23hat = new Vector3D();
        v23hat.E(v23);
        v23hat.normalize();

        Vector m1 = new Vector3D();
        m1.E(n1Cross);
        m1.XE(v23hat);

        double x = n1Cross.dot(n2Cross);
        double y = m1.dot(n2Cross);

        double phi = Math.atan2(y, x);
        System.out.println(phi);
        System.out.println(Degree.UNIT.fromSim(phi));
    }
}

