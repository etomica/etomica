package etomica.GasMOP.move;

import etomica.action.MoleculeActionTranslateTo;
import etomica.box.RandomPositionSource;
import etomica.box.RandomPositionSourceRectangular;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeArrayList;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.util.random.IRandom;

public class MCMoveGrapheneWrinkle extends MCMoveBox {
    protected final PotentialCompute potentialCompute;
    protected final MoleculeArrayList reservoir;
    protected final MoleculeActionTranslateTo atomTranslator;
    //chemical potential
    protected double mu;
    protected ISpecies speciesGraphene;
    protected IMolecule testMolecule;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected boolean insert;
    protected IMoleculeList moleculeList;
    protected IRandom random;
    protected RandomPositionSource positionSource;

    public MCMoveGrapheneWrinkle (PotentialCompute potentialCompute, IRandom random, Space space){
        super();
        this.potentialCompute = potentialCompute;
        atomTranslator = new MoleculeActionTranslateTo(space);
        reservoir = new MoleculeArrayList();
        this.random = random;
        perParticleFrequency = true;
        positionSource = new RandomPositionSourceRectangular(space, random);
    }

    public void setSpecies(ISpecies s) {
        speciesGraphene = s;
        if (box != null) {
            moleculeList = box.getMoleculeList(speciesGraphene);
        }
    }

    public boolean doTrial(){


        return true;
    }

    public double getChi(double temperature){
        uNew = potentialCompute.computeOneOldMolecule(testMolecule);
        return Math.exp(-(uNew - uOld) / temperature);
    }

    public void acceptNotify(){

    }

    public void rejectNotify(){

    }

    public double energyChange() {
        return uNew - uOld;
    }

}
