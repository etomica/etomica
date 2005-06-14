package etomica.integrator.mcmove;

import etomica.Phase;
import etomica.PotentialMaster;
import etomica.action.AtomActionTranslateBy;
import etomica.action.AtomGroupAction;
import etomica.atom.AtomSourceRandomMolecule;
import etomica.space.Vector;

/**
 * Standard Monte Carlo molecule-displacement trial move.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 7/9/02 Added energyChange() method
  */
public class MCMoveMolecule extends MCMoveAtom {
    
    protected final AtomGroupAction moveMoleculeAction;
    protected final Vector groupTranslationVector;

    public MCMoveMolecule(PotentialMaster potentialMaster) {
        super(potentialMaster);
        AtomActionTranslateBy translator = new AtomActionTranslateBy(potentialMaster.getSpace());
        groupTranslationVector = translator.getTranslationVector();
        moveMoleculeAction = new AtomGroupAction(translator);
        
        //set directive to exclude intramolecular contributions to the energy

        //TODO enable meter to do this
        //       iteratorDirective.addCriterion(new IteratorDirective.PotentialCriterion() {
 //           public boolean excludes(Potential p) {return (p instanceof Potential1.Intramolecular);}
 //       });
        setName("MCMoveMolecule");
        setAtomSource(new AtomSourceRandomMolecule());
    }
    

    public boolean doTrial() {
        Phase phase = phases[0];
        if(phase.moleculeCount()==0) return false;
        
        atom = atomSource.getAtom();

        energyMeter.setTarget(atom);
        uOld = energyMeter.getDataAsScalar(phase);
        groupTranslationVector.setRandomCube();
        groupTranslationVector.TE(stepSize);
        moveMoleculeAction.actionPerformed(atom);
        uNew = Double.NaN;
        return true;
    }//end of doTrial
    
    public void rejectNotify() {
        groupTranslationVector.TE(-1);
        moveMoleculeAction.actionPerformed(atom);
    }
        
/*    public static void main(String[] args) {
        
	    IntegratorMC integrator = new IntegratorMC();
	    MCMoveMolecule mcMove = new MCMoveMolecule(integrator);
	    SpeciesSpheres species = new SpeciesSpheres(20,3);
	    Phase phase = new Phase();
	    P2LennardJones potential = new P2LennardJones();
	    Controller controller = new Controller();
	    DisplayPhase display = new DisplayPhase();

		MeterEnergy energy = new MeterEnergy();
		energy.setPhase(phase);
		energy.setHistorying(true);
		energy.setActive(true);		
		energy.getHistory().setNValues(500);		
		DisplayPlot plot = new DisplayPlot();
		plot.setLabel("Energy");
		plot.setDataSources(energy.getHistory());
		
		integrator.setSleepPeriod(2);
		
		Simulation.instance.elementCoordinator.go();
	    
        potential.setIterator(new AtomPairIterator(phase));
        potential.set(species.getAgent(phase));
        
        Simulation.makeAndDisplayFrame(Simulation.instance);
    }//end of main
  */  
}