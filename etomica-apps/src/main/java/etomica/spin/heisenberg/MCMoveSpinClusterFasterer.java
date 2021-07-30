package etomica.spin.heisenberg;

import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.box.Box;
import etomica.integrator.IntegratorMCFasterer;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.NeighborManager;
import etomica.potential.compute.PotentialCompute;
import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

import java.util.HashSet;

/**
 * Cluster move as described by Wolff in PRL v62 361 (1989)
 */
public class MCMoveSpinClusterFasterer extends MCMoveBox {

    protected final Space space;
    protected final PotentialCompute potentialCompute;
    protected final NeighborIterator nbrItertator;
    protected final IRandom random;
    protected final HashSet<IAtom> clusterAtoms, jNeighbors;
    protected final Vector reflectionVector;
    protected final AtomArrayList atomPairs;
    protected final IntegratorMCFasterer integratorMC;
    protected double J;
    protected final AtomIteratorArrayListSimple atomIterator;

    public MCMoveSpinClusterFasterer(Space space, IRandom random, PotentialCompute potentialCompute, NeighborManager nbrManager, IntegratorMCFasterer integratorMC, double J) {
        super(null);
        this.potentialCompute = potentialCompute;
        this.space = space;
        this.random = random;
        clusterAtoms = new HashSet<>();
        jNeighbors = new HashSet<>();
        reflectionVector = space.makeVector();
        atomPairs = new AtomArrayList();
        this.integratorMC = integratorMC;
        nbrItertator = nbrManager.makeNeighborIterator();
        this.J = J;
        atomIterator = new AtomIteratorArrayListSimple();
    }

    @Override
    public AtomIterator affectedAtoms() {
        AtomArrayList atomList = (AtomArrayList) atomIterator.getList();
        atomList.clear();
        atomList.addAll(clusterAtoms);
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
        IAtom atomI = atoms.get(random.nextInt(atoms.size()));
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
            while (atomPairs.size() > 0) {
                IAtom atomJ = atomPairs.remove(atomPairs.size() - 1);
                atomI = atomPairs.remove(atomPairs.size() - 1);
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
        nbrItertator.iterAllNeighbors(atomI.getLeafIndex(), new NeighborIterator.NeighborConsumer() {
            @Override
            public void accept(IAtom jAtom, Vector rij) {
                if (!clusterAtoms.contains(jAtom)) jNeighbors.add(jAtom);
            }
        });
    }

    public void setBox(Box box) {
        super.setBox(box);
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
        potentialCompute.computeAll(false);
    }

    @Override
    public void rejectNotify() {
        throw new RuntimeException("Should never be rejected");
    }
}
