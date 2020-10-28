package etomica.potential.compute;

import etomica.atom.IAtom;
import etomica.molecule.IMolecule;
import etomica.space.Vector;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class PotentialComputeAggregate implements PotentialCompute {
    private final List<PotentialCompute> potentialComputes;

    public PotentialComputeAggregate() {
        potentialComputes = new ArrayList<>(1);
    }

    public PotentialComputeAggregate(PotentialCompute... computes) {
        this();
        this.add(computes);
    }

    public void add(PotentialCompute... compute) {
        this.potentialComputes.addAll(Arrays.asList(compute));
    }

    @Override
    public void init() {
        this.potentialComputes.forEach(PotentialCompute::init);
    }

    @Override
    public Vector[] getForces() {
        return this.potentialComputes.get(0).getForces();
    }

    @Override
    public double getLastVirial() {
        return this.potentialComputes.stream().mapToDouble(PotentialCompute::getLastVirial).sum();
    }

    @Override
    public double getOldEnergy() {
        return this.potentialComputes.stream().mapToDouble(PotentialCompute::getOldEnergy).sum();
    }

    @Override
    public void updateAtom(IAtom atom) {
        for (PotentialCompute c : this.potentialComputes) {
            c.updateAtom(atom);
        }
    }

    @Override
    public double computeAll(boolean doForces) {
        double sum = 0;
        for (PotentialCompute potentialCompute : this.potentialComputes) {
            sum += potentialCompute.computeAll(doForces);
        }
        if (doForces) {
            Vector[] forces = potentialComputes.get(0).getForces();
            potentialComputes.stream()
                    .skip(1)
                    .forEach(compute -> {
                        Vector[] f = compute.getForces();
                        for (int i = 0; i < forces.length; i++) {
                            forces[i].PE(f[i]);
                        }
                    });
        }
        return sum;
    }

    @Override
    public double computeOneOld(IAtom iAtom) {
        return this.potentialComputes.stream().mapToDouble(c -> c.computeOneOld(iAtom)).sum();
    }

    @Override
    public double computeOneOldMolecule(IMolecule molecule) {
        return this.potentialComputes.stream().mapToDouble(c -> c.computeOneOldMolecule(molecule)).sum();
    }

    @Override
    public double computeOne(IAtom iAtom) {
        return this.potentialComputes.stream().mapToDouble(c -> c.computeOne(iAtom)).sum();
    }

    @Override
    public double computeOneMolecule(IMolecule molecule) {
        return this.potentialComputes.stream().mapToDouble(c -> c.computeOneMolecule(molecule)).sum();
    }

    @Override
    public void processAtomU(double fac) {
        for (PotentialCompute compute : this.potentialComputes) {
            compute.processAtomU(fac);
        }
    }
}
