package etomica.integrator.mcmove;
import etomica.action.AtomTransform;
import etomica.atom.Atom;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorNull;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.atom.iterator.AtomIteratorTree;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.RotationTensor;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Monte Carlo trial that rotates the atoms of a molecule about its first atom.
 *
 * Has a bug, probably associated with incorrect replacement of the molecule when
 * rejecting the trial.
 */
 
 /* History of changes
  * 7/9/02 Added energyChange() method
  * 10/06/03 (DAK) added check in constructor to ensure simulation is 2D
  */
  
public class MCMoveRotateMolecule extends MCMoveStep {
    
    private final MeterPotentialEnergy energyMeter;
    private final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    
    //not sure that tree iterator will give all leaf atoms of a molecule
    //because it assumes all subgroups have same structure
    private final AtomIterator leafAtomIterator = new AtomIteratorTree();
    
    private transient double uOld;
    private transient double uNew = Double.NaN;
    private transient Atom molecule;
    private transient Vector r0;
    private transient RotationTensor rotationTensor;

    public MCMoveRotateMolecule(PotentialMaster potentialMaster, Space space) {
        super(potentialMaster, new MCMoveStepTracker(), 1);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        if(space.D() != 2) throw new RuntimeException("MCMoveRotateMolecule suitable only for 2-D simulation");
        rotationTensor = space.makeRotationTensor();
        r0 = space.makeVector();
        setStepSizeMax(Math.PI);
        setStepSizeMin(0.0);
        setStepSize(Math.PI/2.0);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);
        //set directive to exclude intramolecular contributions to the energy
        //TODO enable this
//        iteratorDirective.addCriterion(new IteratorDirective.PotentialCriterion() {
//            public boolean excludes(Potential p) {return (p instanceof Potential1.Intramolecular);}
//        });
    }
     
    public boolean doTrial() {
        if(phase.moleculeCount()==0) {molecule = null; return false;}
        molecule = phase.randomMolecule();
        energyMeter.setTarget(molecule);
        uOld = energyMeter.getDataAsScalar();
        //update for 3D
        double dTheta = (2*Simulation.random.nextDouble() - 1.0)*stepSize;
        rotationTensor.setAxial(2,dTheta);
       // molecule.coord.transform(molecule.coord.position(), rotationTensor);
        leafAtomIterator.reset();
        r0.E(molecule.node.firstLeafAtom().coord.position());
//        AtomTransform.doAction(leafAtomIterator, molecule.coord.position(), rotationTensor);
        AtomTransform.doTransform(leafAtomIterator, r0, rotationTensor);
        uNew = Double.NaN;
        return true;
    }//end of doTrial
    
    public double getA() {return 1.0;}
    
    public double getB() {
        energyMeter.setTarget(molecule);
        uNew = energyMeter.getDataAsScalar();
        return -(uNew - uOld);
    }
    
    public void acceptNotify() {  /* do nothing */}
    
    public void rejectNotify() {
        rotationTensor.invert();
        leafAtomIterator.reset();
        AtomTransform.doTransform(leafAtomIterator, r0, rotationTensor);
    }
 
    public double energyChange(Phase p) {return (p == phase) ? uNew - uOld : 0.0;}
    
    public final AtomIterator affectedAtoms(Phase p) {
        if(p != phase) return AtomIteratorNull.INSTANCE;
        affectedAtomIterator.setAtom(molecule);
        affectedAtomIterator.reset();
        return affectedAtomIterator;
    }

    public void setPhase(Phase p) {
        super.setPhase(p);
        energyMeter.setPhase(p);
    }
}
