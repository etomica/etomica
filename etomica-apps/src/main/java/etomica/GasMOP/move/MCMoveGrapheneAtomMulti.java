package etomica.GasMOP.move;

import etomica.atom.AtomSource;
import etomica.atom.AtomSourceRandomLeaf;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.molecule.IMolecule;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Vector;
import etomica.util.random.IRandom;
import etomica.virial.BoxCluster;

import java.util.ArrayList;
import java.util.List;

public class MCMoveGrapheneAtomMulti  extends MCMoveBoxStep {
    protected final IRandom random;
    protected final PotentialCompute potentialCompute;
    protected final List<Vector> oldPositions = new ArrayList<>();
    protected IMolecule molecule;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected Vector[] translationVectors;
    protected int maxAtoms;

    protected AtomSource[] atomSourceArray;
    protected int startAtom = 1;
    protected boolean imposePBC = false;

    public MCMoveGrapheneAtomMulti(IRandom random, PotentialCompute potentialCompute, Box box, int maxAtoms) {
        super();
        this.random = random;
        this.potentialCompute = potentialCompute;
        setStepSize(1.2);
        setBox(box);
    }

    public void setBox(Box p) {
        super.setBox(p);
        if (translationVectors == null) {
            translationVectors = new Vector[box.getLeafList().size() - startAtom];
            for (int i=0; i<translationVectors.length; i++) {
                translationVectors[i] = box.getSpace().makeVector();
            }
        }
    }

    public void setStartAtom(int newStartAtom) {
        startAtom = newStartAtom;
        if (translationVectors != null && translationVectors.length != box.getLeafList().size() - startAtom) {
            translationVectors = new Vector[box.getLeafList().size() - startAtom];
            for (int i = 0; i < translationVectors.length; i++) {
                translationVectors[i] = box.getSpace().makeVector();
            }
        }

    }
    public boolean doTrial() {
        IAtomList leafAtoms = box.getLeafList();
        uOld = potentialCompute.computeOneOldMolecule(molecule);
        for(int i = startAtom; i<leafAtoms.size(); i++) {
            translationVectors[i - startAtom].setRandomCube(random);
            translationVectors[i - startAtom].TE(stepSize);
            Vector r = leafAtoms.get(i).getPosition();
            r.PE(translationVectors[i - startAtom]);
            if (imposePBC) r.PE(box.getBoundary().centralImage(r));
        }
        uNew = potentialCompute.computeOneMolecule(molecule);
        ((BoxCluster)box).trialNotify();
        return true;
    }

    public int getStartAtom() {
        return startAtom;
    }

    public void setDoImposePBC(boolean doImposePBC) {
        imposePBC = doImposePBC;
    }

    public double getChi(double temperature) {
        uNew = potentialCompute.computeOneMolecule(molecule);
        return Math.exp(-(uNew - uOld) / temperature);
    }

    public double energyChange() {
        return uNew - uOld;
    }

    public void rejectNotify() {
        IAtomList leafAtoms = box.getLeafList();
        for(int i = startAtom; i<leafAtoms.size(); i++) {
            Vector r = leafAtoms.get(i).getPosition();
            r.ME(translationVectors[i - startAtom]);
            if (imposePBC) r.PE(box.getBoundary().centralImage(r));
        }
        ((BoxCluster)box).rejectNotify();

        molecule.getChildList().forEach(atom -> {
            atom.getPosition().E(oldPositions.get(atom.getIndex()));
            potentialCompute.updateAtom(atom);
        });
    }

    public void acceptNotify() {
        ((BoxCluster)box).acceptNotify();
    }
}
