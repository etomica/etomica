package simulate;

import java.util.Random;

public class MCMoveInsertDelete extends MCMove {
    
    private transient final double[] dr = new double[Space.D];
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
    
    private final void trialInsert(Phase phase) {
        Species s = phase.firstSpecies();
        double uNew;
        phase.space.randomVector(dr, rand);  //random point in volume
        Molecule m = new Molecule(s,s.nAtomsPerMolecule);  //makes molecule without adding to species
        m.translate(dr);
        uNew = phase.potentialEnergy.insertionValue(m);
        if(uNew == Double.MAX_VALUE) {return;}        //overlap
        double bNew = Math.exp((mu-uNew)/parentIntegrator.temperature)*phase.space.volume/(s.nMolecules+1);
        if(bNew > 1.0 || bNew > rand.nextDouble()) {  //accept
            s.addMolecule(m);
        }           
    }
    
    private final void trialDelete(Phase phase) {
        Species s = phase.firstSpecies();
        if(s.nMolecules == 0) {return;}
        double bOld, bNew;
        int i = (int)(rand.nextDouble()*s.nMolecules);
        Molecule m = s.firstMolecule;
        for(int j=i; --j>=0; ) {m = m.getNextMolecule();}
        bOld = Math.exp((mu-phase.potentialEnergy.currentValue(m))/parentIntegrator.temperature);
        bNew = s.nMolecules/(phase.space.volume);
        if(bNew > bOld || bNew > rand.nextDouble()*bOld) {  //accept
            s.deleteMolecule(m);
        }           
    }

    public final void setMu(double mu) {this.mu = mu;}
    public final double getMu() {return mu;}
}