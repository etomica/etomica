package etomica.spin.heisenberg;

import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.box.Box;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.nbr.site.Api1ASite;
import etomica.nbr.site.PotentialMasterSite;
import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

import java.util.HashSet;

/**
 * Cluster move as described by Wolff in PRL v62 361 (1989)
 */
public class MCMoveSpinCluster extends MCMoveBox {

    protected final Space space;
    protected final IRandom random;
    protected final HashSet<IAtom> clusterAtoms, jNeighbors;
    protected final Vector reflectionVector;
    protected final AtomArrayList atomPairs;
    protected final IntegratorMC integratorMC;
    protected final Api1ASite neighbors;
    protected double J;
    protected final AtomIteratorArrayListSimple atomIterator;

    public MCMoveSpinCluster(Space space, IRandom random, PotentialMasterSite potentialMaster, IntegratorMC integratorMC, double J) {
        super(potentialMaster);
        this.space = space;
        this.random = random;
        clusterAtoms = new HashSet<>();
        jNeighbors = new HashSet<>();
        reflectionVector = space.makeVector();
        atomPairs = new AtomArrayList();
        this.integratorMC = integratorMC;
        neighbors = new Api1ASite(space.D(), potentialMaster.getCellAgentManager());
        neighbors.setDirection(null);
        this.J = J;
        atomIterator = new AtomIteratorArrayListSimple();
    }

    @Override
    public AtomIterator affectedAtoms() {
        AtomArrayList atomList = (AtomArrayList) atomIterator.getList();
        atomList.clear();
        for (IAtom atom : clusterAtoms) {
            atomList.add(atom);
        }
        return atomIterator;
    }

    @Override
    public double energyChange() {
        // wrong
        return 0;
    }

    @Override
    public boolean doTrial() {
        double temperature = integratorMC.getTemperature();
        reflectionVector.setRandomSphere(random);
        IAtomList atoms = box.getLeafList();
        IAtom atomI = atoms.getAtom(random.nextInt(atoms.getAtomCount()));
        clusterAtoms.clear();
        clusterAtoms.add(atomI);
        atomPairs.clear();
        outer:
        while (true) {
            jNeighbors.clear();
            gatherNeighbors(atomI);
            rotateAtom(atomI);
            for (IAtom atomJ : jNeighbors) {
                atomPairs.add(atomI);
                atomPairs.add(atomJ);
            }
            while (atomPairs.getAtomCount() > 0) {
                IAtom atomJ = atomPairs.remove(atomPairs.getAtomCount() - 1);
                atomI = atomPairs.remove(atomPairs.getAtomCount() - 1);
                if (clusterAtoms.contains(atomJ)) continue;
                Vector oi = ((IAtomOriented) atomI).getOrientation().getDirection();
                Vector oj = ((IAtomOriented) atomJ).getOrientation().getDirection();
                double x = reflectionVector.dot(oi) * reflectionVector.dot(oj);
                if (x > 0) continue;
                double p = 1 - Math.exp(2 * J * x / temperature);
                if (p < random.nextDouble()) continue;
                atomI = atomJ;
                clusterAtoms.add(atomI);
                continue outer;
            }
            break;
        }
        return true;
    }

    private void gatherNeighbors(IAtom atomI) {
        neighbors.setTarget(atomI);
        neighbors.reset();
        for (IAtomList pair = neighbors.next(); pair != null; pair = neighbors.next()) {
            IAtom atomJ = pair.getAtom(0);
            if (atomJ == atomI) atomJ = pair.getAtom(1);
            if (clusterAtoms.contains(atomJ)) continue;
            jNeighbors.add(atomJ);
        }
    }

    public void setBox(Box box) {
        super.setBox(box);
        neighbors.setBox(box);
    }

    private void rotateAtom(IAtom atom) {
        IOrientation o = ((IAtomOriented) atom).getOrientation();
        Vector direction = o.getDirection();
        direction.PEa1Tv1(-2 * direction.dot(reflectionVector), reflectionVector);
        o.setDirection(direction);
    }

    @Override
    public double getChi(double temperature) {
        return 1;
    }

    @Override
    public void acceptNotify() {
    }

    @Override
    public void rejectNotify() {
        throw new RuntimeException("Should never be rejected");
    }
}
