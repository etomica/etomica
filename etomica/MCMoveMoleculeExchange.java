package etomica;

import java.util.Random;

/**
 * Performs a trial that results in the exchange of a molecule from one phase to another.
 * Primary use is as an elementary move in a Gibbs ensemble simulation
 *
 * @author David Kofke
 */
public final class MCMoveMoleculeExchange extends MCMove {
    
    private final Random rand = new Random();
    private Phase firstPhase;
    private Phase secondPhase;
    private final double ROOT;

    public MCMoveMoleculeExchange(Integrator parent) {
        super();
        parentIntegrator = parent;
        ROOT = 1.0/(double)parentIntegrator.parentSimulation().space().D();
        setTunable(false);
        setPerParticleFrequency(true);
    }
    
    /**
     * Overrides superclass method so that it performs no action.
     * Must set using method that takes an array of phases.
     */
    public void setPhase(Phase p) {}

    public void setPhase(Phase[] p) {
        firstPhase = p[0];
        secondPhase = p[1];
    }
    
    /**
     * Overrides superclass method so that frequency is set based on number of molecules in both phases.
     */
    public void resetFrequency() {
        frequency = perParticleFrequency ? nominalFrequency*(firstPhase.moleculeCount+secondPhase.moleculeCount) : nominalFrequency;
    }
    
    //under revision--- does not work for multiatomics, since intramolecular energy is not considered
    public void thisTrial() {
        Phase iPhase, dPhase;
        if(java.lang.Math.random() < 0.5) {
            iPhase = firstPhase;
            dPhase = secondPhase;
        }
        else {
            iPhase = secondPhase;
            dPhase = firstPhase;
        }
        if(dPhase.moleculeCount() == 0) {return;} //no molecules to delete; trial is over
        
        Molecule m = dPhase.randomMolecule();  //select random molecule to delete
        Species species = m.parentSpecies();
        
        Species.Agent iSpecies = species.getAgent(iPhase);  //insertion-phase speciesAgent
        Species.Agent dSpecies = species.getAgent(dPhase);  //deletion-phase species Agent
        
        double uOld = dPhase.energy.meterPotential().currentValue(m); //get its contribution to energy

        iPhase.addMolecule(m,iSpecies);
        m.displaceTo(iPhase.randomPosition());         //place at random in insertion phase
        double uNew = iPhase.energy.meterPotential().currentValue(m); //get its new energy
        if(uNew == Double.MAX_VALUE) {  //overlap; reject
            m.replace();                //put it back 
            dPhase.addMolecule(m,dSpecies);
            return;        
        }
        
        double bFactor = (dSpecies.nMolecules+1)/dPhase.volume()  //acceptance probability
                         * iPhase.volume()/(iSpecies.nMolecules)                   //note that dSpecies.nMolecules has been decremented
                         * Math.exp(-(uNew-uOld)/parentIntegrator.temperature);    //and iSpecies.nMolecules has been incremented
        if(bFactor > 1.0 || bFactor > rand.nextDouble()) {  //accept
            nAccept++;
        }
        else {              //reject
            m.replace();    //put it back
            dPhase.addMolecule(m,dSpecies); //place molecule in phase
        }
    }//end of thisTrial
}//end of MCMoveMoleculeExchange