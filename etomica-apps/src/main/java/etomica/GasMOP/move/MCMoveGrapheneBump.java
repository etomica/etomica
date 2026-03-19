package etomica.GasMOP.move;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculeSource;
import etomica.molecule.MoleculeSourceRandomMolecule;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import etomica.util.random.IRandom;

import java.util.ArrayList;
import java.util.List;
/* Δz(x,y) = A exp ( ((x−x0)^2+(y−y0)^2) / R ^ 2)*/

public class MCMoveGrapheneBump  extends MCMoveBoxStep {
    protected final PotentialCompute potentialCompute;
    protected final List<Vector> oldPositions = new ArrayList<>();
    protected final IRandom random;
    protected IMolecule molecule;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected MoleculeSource moleculeSource;
    protected Space space;
    protected double A, R;
    protected double Aval, Rval, dz, val, xn, yn, x0, y0;
    public MCMoveGrapheneBump(IRandom random, PotentialCompute potentialCompute, Box box){
        super();
        this.potentialCompute = potentialCompute;
        this.random = random;
        this.space = box.getSpace();
        MoleculeSourceRandomMolecule source = new MoleculeSourceRandomMolecule();
        source.setRandomNumberGenerator(random);
        moleculeSource = source;
        setMaxBumpSize(3, 3);
        setBox(box);
        setBoxSize(box);
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
        Aval = A * random.nextDouble();
        Rval = R * random.nextDouble();
        return true;
    }

    protected void doTransform(){
        IAtomList childList = molecule.getChildList();
        //Δz(x,y) = A exp ( ((x−x0)^2+(y−y0)^2) / R ^ 2)
        for (IAtom a : childList) {
            Vector r = a.getPosition();
            xn = (r.getX(0) - x0);
            yn = (r.getX(1) - y0);
            val = (xn * xn + yn * yn) / (Rval * Rval);
            dz = Aval * Math.exp(val);
            r.PE(new Vector3D(0, 0, dz));
            r.PE(box.getBoundary().centralImage(r));
            potentialCompute.updateAtom(a);
        }
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

    public void setBoxSize(Box box){
        double dimX = box.getBoundary().getBoxSize().getX(0);
        double dimY = box.getBoundary().getBoxSize().getX(0);
        this.x0 = 0.6 * dimX;
        this.y0 = 0.6 * dimY;
    }

    public void setMoleculeSource(MoleculeSource source) {
        moleculeSource = source;
    }

    public void setMaxBumpSize(double maxBumpHeight, double maxBumpSize){this.A = maxBumpHeight; this.R = maxBumpSize;}
}

