package etomica;

import java.util.Random;

public class MCMoveInsertDelete extends MCMove {
    
    private final Random rand = new Random();
    double mu;
    private final IteratorDirective iteratorDirective = new IteratorDirective();
    private final PotentialCalculation.EnergySum energy = new PotentialCalculation.EnergySum();

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
        energy.reset();
        m.atomIterator.reset();
        while(m.atomIterator.hasNext()) {
            phase.potential().calculate(iteratorDirective, energy);
        }
        uNew = energy.sum();
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
        energy.reset();
        m.atomIterator.reset();
        while(m.atomIterator.hasNext()) {
            phase.potential().calculate(iteratorDirective, energy);
        }
        double uOld = energy.sum();
        bOld = Math.exp((mu-uOld)/parentIntegrator.temperature);
        bNew = s.nMolecules/(phase.volume());
        if(bNew > bOld || bNew > rand.nextDouble()*bOld) {  //accept
            phase.deleteMolecule(m);
            nAccept++;
        }           
    }

    public final void setMu(double mu) {this.mu = mu;}
    public final double getMu() {return mu;}
    public final etomica.units.Dimension getMuDimension() {return etomica.units.Dimension.ENERGY;}
    
    public static void main(String[] args) {
        etomica.simulations.HsMc2d sim = new etomica.simulations.HsMc2d();
        Simulation.instance = sim;

        MeterNMolecules meterN = new MeterNMolecules();
        DisplayBox box = new DisplayBox(meterN);
        box.setUpdateInterval(10);

		Simulation.instance.elementCoordinator.go();

        MCMoveInsertDelete mcMoveInsDel = new MCMoveInsertDelete();
        sim.integrator.add(mcMoveInsDel);
		                                    
        Simulation.makeAndDisplayFrame(sim);
    }//end of main
}//end of MCMoveInsertDelete