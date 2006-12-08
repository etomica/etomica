package etomica.integrator.mcmove;

import etomica.action.AtomActionTranslateBy;
import etomica.action.AtomGroupAction;
import etomica.atom.AtomSourceRandomMolecule;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Vector;

/**
 * Standard Monte Carlo molecule-displacement trial move.
 *
 * @author David Kofke
 */
public class MCMoveMolecule extends MCMoveAtom {
    
    private static final long serialVersionUID = 1L;
    protected final AtomGroupAction moveMoleculeAction;
    protected final Vector groupTranslationVector;

    public MCMoveMolecule(Simulation sim) {
        this(sim.getPotentialMaster(),sim.getDefaults().atomSize,sim.getDefaults().boxSize*0.5,
                sim.getDefaults().ignoreOverlap);
    }
    
    public MCMoveMolecule(PotentialMaster potentialMaster, double stepSize,
            double stepSizeMax, boolean ignoreOverlap) {
        super(potentialMaster,stepSize,stepSizeMax,ignoreOverlap);
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
        if(phase.moleculeCount()==0) return false;
        
        atom = atomSource.getAtom();

        energyMeter.setTarget(atom);
        uOld = energyMeter.getDataAsScalar();
        groupTranslationVector.setRandomCube();
        groupTranslationVector.TE(stepSize);
        moveMoleculeAction.actionPerformed(atom);
        uNew = Double.NaN;
        return true;
    }
    
    public void rejectNotify() {
        groupTranslationVector.TE(-1);
        moveMoleculeAction.actionPerformed(atom);
    }
        
}