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
        setTunable(false);
    }
    
    public final void thisTrial() {
        if(rand.nextDouble() < 0.5) {
            trialInsert();  //update for more species
        }
        else {
            trialDelete();
        }
    }
                                                                                                                                                                                                                                                                                                                                                                       
    //Not suited to mixtures
    private final void trialInsert() {
        Species.Agent s = phase.firstSpecies();
        double uNew;
        Molecule m = s.parentSpecies().getMolecule();  //makes molecule without adding to species
        phase.addMolecule(m,s);
        m.translateTo(phase.randomPosition());
        uNew = phase.energy.meterPotential().currentValue(m);
        if(uNew == Double.MAX_VALUE) {  //overlap
            phase.deleteMolecule(m);
            return;
        }      
        double bNew = Math.exp((mu-uNew)/parentIntegrator.temperature)*phase.volume()/(s.moleculeCount()+1);
        if(bNew < 1.0 && bNew < rand.nextDouble()) {  //reject
            phase.deleteMolecule(m);
        }
        else nAccept++;
    }
    
    //Not suited to mixtures
    private final void trialDelete() {
        Species.Agent s = phase.firstSpecies();
        if(s.moleculeCount() == 0) {return;}
        double bOld, bNew;
        int i = (int)(rand.nextDouble()*s.moleculeCount());
        Molecule m = s.firstMolecule;
        for(int j=i; --j>=0; ) {m = m.nextMolecule();}
        bOld = Math.exp((mu-phase.energy.meterPotential().currentValue(m))/parentIntegrator.temperature);
        bNew = s.nMolecules/(phase.volume());
        if(bNew > bOld || bNew > rand.nextDouble()*bOld) {  //accept
            phase.deleteMolecule(m);
            nAccept++;
        }           
    }

    public final void setMu(double mu) {this.mu = mu;}
    public final double getMu() {return mu;}
}