package etomica.integrator.mcmove;

import etomica.Atom;
import etomica.AtomFactory;
import etomica.AtomIterator;
import etomica.AtomTreeNodeGroup;
import etomica.Phase;
import etomica.PotentialMaster;
import etomica.Simulation;
import etomica.Species;
import etomica.SpeciesAgent;
import etomica.action.AtomActionTranslateTo;
import etomica.atom.AtomList;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.MCMove;

/**
 * Elementary Monte Carlo move in which a molecule of a specified species is
 * inserted into or removed from a phase.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 07/09/02 (DAK) Added energyChange() method
  * 09/19/02 (DAK) Minor change in doTrial for case were deleting with N = 0
  */
public class MCMoveInsertDelete extends MCMove {
    
    //chemical potential
    protected double mu;
    
    //directive must specify "BOTH" to get energy with all atom pairs
    protected final MeterPotentialEnergy energyMeter;
	protected Species species;
	protected AtomTreeNodeGroup speciesAgentNode;
    protected SpeciesAgent speciesAgent;
	protected final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
	protected Atom testMolecule;
	protected double uOld;
	protected double uNew = Double.NaN;
	protected boolean insert;
	protected final AtomList reservoir;
    protected final AtomActionTranslateTo atomTranslator;
    protected AtomFactory moleculeFactory;

    public MCMoveInsertDelete(PotentialMaster potentialMaster) {
        super(potentialMaster, 1);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.10);
        setMu(0.0);
        setTunable(false);
        energyMeter.setIncludeLrc(true);
        atomTranslator = new AtomActionTranslateTo(potentialMaster.getSpace());
        reservoir = new AtomList();
    }
    
//perhaps should have a way to ensure that two instances of this class aren't assigned the same species
    public void setSpecies(Species s) {
        species = s;
        if(phases[0] != null) {
            speciesAgent = species.getAgent(phases[0]);
            speciesAgentNode = (AtomTreeNodeGroup)speciesAgent.node; 
        }
        moleculeFactory = species.moleculeFactory();
    }
    public Species getSpecies() {return species;}
    
    public void setPhase(Phase[] p) {
        super.setPhase(p);
        if(species != null) {
            speciesAgent = species.getAgent(phases[0]);
            speciesAgentNode = (AtomTreeNodeGroup)speciesAgent.node; 
        }
    }
    
    /**
     * Chooses and performs with equal probability an elementary molecule insertion
     * or deletion.
     */
    public boolean doTrial() {
        insert = Simulation.random.nextDouble() < 0.5;
        if(insert) {
            uOld = 0.0;
            
            if(!reservoir.isEmpty()) testMolecule = reservoir.removeFirst();
            else testMolecule = moleculeFactory.makeAtom();
            phases[0].addMolecule(testMolecule, speciesAgent);

            atomTranslator.setDestination(phases[0].randomPosition());
            atomTranslator.actionPerformed(testMolecule);
        } else {//delete
            if(speciesAgent.moleculeCount() == 0) {
                testMolecule = null;
                return false;
            }
            testMolecule = speciesAgent.randomMolecule();
            //delete molecule only upon accepting trial
            energyMeter.setTarget(testMolecule);
            uOld = energyMeter.getDataAsScalar(phases[0]);
        } 
        uNew = Double.NaN;
        return true;
    }//end of doTrial
    
    public double lnTrialRatio() {//note that moleculeCount() gives the number of molecules after the trial is attempted
        return insert ? Math.log(phases[0].volume()/speciesAgent.moleculeCount()) 
                      : Math.log((speciesAgent.moleculeCount()+1)/phases[0].volume());        
    }
    
    public double lnProbabilityRatio() {
        if(insert) {
            energyMeter.setTarget(testMolecule);
            uNew = energyMeter.getDataAsScalar(phases[0]);
           return (+mu - uNew)/temperature;
        }
        uNew = 0.0;
        return (-mu + uOld)/temperature;
    }
    
    public void acceptNotify() {
        //      accepted deletion - remove from phase and add to reservoir 
        if(!insert) {
            phases[0].removeMolecule(testMolecule);
            reservoir.add(testMolecule.seq);
        }
    }
    
    public void rejectNotify() {
        //      rejected insertion - remove from phase and return to reservoir
        if(insert) {
            phases[0].removeMolecule(testMolecule);
            reservoir.add(testMolecule.seq);
        }
    }
    
    public double energyChange(Phase phase) {return (this.phases[0] == phase) ? uNew - uOld : 0.0;}

    /**
     * Returns an iterator giving molecule that is being added or deleted 
     * in the current or most recent trial.
     */
    public final AtomIterator affectedAtoms(Phase phase) {
        if(this.phases[0] != phase) return AtomIterator.NULL;
        affectedAtomIterator.setAtom(testMolecule);
        return affectedAtomIterator;
    }

    /**
     * Mutator method for the chemical potential of the insertion/deletion species.
     */
    public final void setMu(double mu) {this.mu = mu;}
    /**
     * Accessor method for the chemical potential of th insertion/deletion species.
     */
    public final double getMu() {return mu;}
    /**
     * Indicates that chemical potential has dimensions of energy.
     */
    public final etomica.units.Dimension getMuDimension() {return etomica.units.Dimension.ENERGY;}
///*    
//    public static void main(String[] args) {
//        Default.TRUNCATE_POTENTIALS = false;
////        etomica.simulations.HsMc2d sim = new etomica.simulations.HsMc2d();
////		MeterNMolecules meterN = new MeterNMolecules();
////		etomica.graphics.DisplayBox box = new etomica.graphics.DisplayBox((DatumSource)meterN);
////		box.setUpdateInterval(10);
////        
////		MCMoveInsertDelete mcMoveInsDel = new MCMoveInsertDelete(sim.integrator);
////		mcMoveInsDel.setSpecies(sim.species);
////		mcMoveInsDel.setMu(-2000.);
////        
////		sim.integrator(0).setTemperature(1.0);
////		                                    
////		etomica.graphics.DeviceSlider slider = new etomica.graphics.DeviceSlider(mcMoveInsDel,"mu");
////		slider.setMinimum(-20);
////		slider.setMaximum(+20);
////		Simulation.instance.elementCoordinator.go();
//		Simulation sim = new etomica.graphics.SimulationGraphic();
//		Controller controller = new Controller();
//
//		SpeciesSpheresMono species = new SpeciesSpheresMono();
//		
//		Phase phase1 = new Phase();
//		Phase phase2 = new Phase();
//		IntegratorMC integrator1 = new IntegratorMC();
//		IntegratorMC integrator2 = new IntegratorMC();
//		phase1.setIntegrator(integrator1);
//		phase2.setIntegrator(integrator2);
//		integrator1.setTemperature(1.0);
//		integrator2.setTemperature(1.0);
//		
//		P2HardSphere potential = new P2HardSphere();
//		potential.setSpecies(species);
//		
//		sim.elementCoordinator.go();
//		
//		MCMoveAtom moveAtom1 = new MCMoveAtom(integrator1);
//		MCMoveAtom moveAtom2 = new MCMoveAtom(integrator2);
//		MCMoveInsertDelete mcMoveInsDel1 = new MCMoveInsertDelete(integrator1);
//		MCMoveInsertDelete mcMoveInsDel2 = new MCMoveInsertDelete(integrator2);
//		mcMoveInsDel1.setSpecies(species);
//		mcMoveInsDel2.setSpecies(species);		
//		mcMoveInsDel1.setMu(-2000.);
//		mcMoveInsDel2.setMu(0.);
//		
//		//TODO this should be one meter and two phases and one meter for each phase
//		MeterNMolecules meterN1 = new MeterNMolecules();
//		meterN1.setPhase(new Phase[] {phase1});
//		etomica.graphics.DisplayBox box1 = new etomica.graphics.DisplayBox((DatumSource)meterN1);
//		box1.setUpdateInterval(10);
//		MeterNMolecules meterN2 = new MeterNMolecules();
//		meterN2.setPhase(new Phase[] {phase2});
//		etomica.graphics.DisplayBox box2 = new etomica.graphics.DisplayBox((DatumSource)meterN2);
//		box1.setUpdateInterval(10);
//        
//		etomica.graphics.DeviceSlider slider = new etomica.graphics.DeviceSlider(mcMoveInsDel1,"mu");
//		slider.setMinimum(-20);
//		slider.setMaximum(+20);
//		
//		etomica.graphics.DisplayPhase display1 = new etomica.graphics.DisplayPhase();
//		etomica.graphics.DisplayPhase display2 = new etomica.graphics.DisplayPhase();
//		display1.setPhase(phase1);
//		display2.setPhase(phase2);
//		Simulation.instance.elementCoordinator.go();
//
//        etomica.graphics.SimulationGraphic.makeAndDisplayFrame(sim);
//    }//end of main
// */   
}//end of MCMoveInsertDelete