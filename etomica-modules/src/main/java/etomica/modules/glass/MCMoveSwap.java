package etomica.modules.glass;

import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIterator;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.util.random.IRandom;

/**
 * MC Move class that swaps the position of two atoms of different species.
 */
public class MCMoveSwap extends MCMoveBox {

    protected IAtom atomA, atomB;
    protected double uOld, uNew;
    protected final ISpecies speciesA, speciesB;
    protected final IRandom random;
    protected MeterPotentialEnergy meterPE;
    protected final Vector tmp;

    public MCMoveSwap(Space space, IRandom random, PotentialMaster potentialMaster, ISpecies speciesA, ISpecies speciesB) {
        super(potentialMaster);
        this.random = random;
        this.speciesA = speciesA;
        this.speciesB = speciesB;
        meterPE = new MeterPotentialEnergy(potentialMaster);
        tmp = space.makeVector();
    }

    @Override
    public void setBox(Box box) {
        super.setBox(box);
        meterPE.setBox(box);
    }

    @Override
    public AtomIterator affectedAtoms() {
        return null;
    }

    @Override
    public double energyChange() {
        return uNew - uOld;
    }

    @Override
    public boolean doTrial() {
        IMoleculeList listA = box.getMoleculeList(speciesA);
        IMoleculeList listB = box.getMoleculeList(speciesB);
        if (listA.size() * listB.size() == 0) return false;
        IMolecule mol = listA.get(random.nextInt(listA.size()));
        atomA = mol.getChildList().get(0);
        mol = listB.get(random.nextInt(listB.size()));
        atomB = mol.getChildList().get(0);

        meterPE.setTarget(atomA);
        uOld = meterPE.getDataAsScalar();
        meterPE.setTarget(atomB);
        uOld += meterPE.getDataAsScalar();
        // since we are swapping, we don't have to worry about the direct
        // interaction being double-counted.  it's unchanged in the new
        // configuration and also double-counted there.
        if (uOld > 1e8) {
            throw new RuntimeException("atom " + atomA + " or " + atomB + " in box " + box + " has an overlap");
        }
        tmp.E(atomA.getPosition());
        atomA.getPosition().E(atomB.getPosition());
        atomB.getPosition().E(tmp);
        if (potential instanceof PotentialMasterList) {
            // need to remove/add so neighbor list manager notices
            box.removeMolecule(atomA.getParentGroup());
            box.addMolecule(atomA.getParentGroup());
            box.removeMolecule(atomB.getParentGroup());
            box.addMolecule(atomB.getParentGroup());
        }
        return true;
    }

    @Override
    public double getChi(double temperature) {
        meterPE.setTarget(atomA);
        uNew = meterPE.getDataAsScalar();
        meterPE.setTarget(atomB);
        uNew += meterPE.getDataAsScalar();
        return Math.exp(-(uNew - uOld) / temperature);
    }

    @Override
    public void acceptNotify() {
//        System.out.println("accepted "+moveTracker.acceptanceProbability());
    }

    @Override
    public void rejectNotify() {
        tmp.E(atomA.getPosition());
        atomA.getPosition().E(atomB.getPosition());
        atomB.getPosition().E(tmp);
        if (potential instanceof PotentialMasterList) {
            box.removeMolecule(atomA.getParentGroup());
            box.addMolecule(atomA.getParentGroup());
            box.removeMolecule(atomB.getParentGroup());
            box.addMolecule(atomB.getParentGroup());
        }
    }
}
