package simulate;

import java.util.Random;

public class MCMoveInsertDelete extends MCMove {
    
    private final Random rand = new Random();
    double mu;

    public MCMoveInsertDelete() {
        super();
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.10);
        setMu(0.0);
    }
    
    public final void thisTrial(Phase phase) {
        if(rand.nextDouble() < 0.5) {
            trialInsert(phase);  //update for more species
        }
        else {
            trialDelete(phase);
        }
    }
    
    //Not suited to mixtures
    private final void trialInsert(Phase phase) {
        Species.Agent s = phase.firstSpecies();
        double uNew;
        Molecule m = s.parentSpecies().makeMolecule();  //makes molecule without adding to species
        m.coordinate.translateToRandom(phase);
        uNew = phase.potentialEnergy.insertionValue(m);
        if(uNew == Double.MAX_VALUE) {return;}        //overlap
        double bNew = Math.exp((mu-uNew)/parentIntegrator.temperature)*phase.volume()/(s.getNMolecules()+1);
        if(bNew > 1.0 || bNew > rand.nextDouble()) {  //accept
            s.addMolecule(m);
        }           
    }
    
    //Not suited to mixtures
    private final void trialDelete(Phase phase) {
        Species.Agent s = phase.firstSpecies();
        if(s.getNMolecules() == 0) {return;}
        double bOld, bNew;
        int i = (int)(rand.nextDouble()*s.getNMolecules());
        Molecule m = s.firstMolecule;
        for(int j=i; --j>=0; ) {m = m.nextMolecule();}
        bOld = Math.exp((mu-phase.potentialEnergy.currentValue(m))/parentIntegrator.temperature);
        bNew = s.nMolecules/(phase.volume());
        if(bNew > bOld || bNew > rand.nextDouble()*bOld) {  //accept
            s.deleteMolecule(m);
        }           
    }

    public final void setMu(double mu) {this.mu = mu;}
    public final double getMu() {return mu;}
}