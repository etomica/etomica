package etomica;
import etomica.action.AtomActionTransform;


public class MCMoveRotateMolecule3D extends MCMove {
    
    private final MeterPotentialEnergy energyMeter;
    private final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    
    private final AtomIterator leafAtomIterator = new AtomIteratorTree();
    
    private transient double uOld;
    private transient double uNew = Double.NaN;
    private transient Atom molecule;
    private transient Space3D.Vector r0;
    private transient Space3D.RotationTensor rotationTensor;
    private int number = 0;
    public int count;
    public int count1;
    public boolean flag = false;
    public boolean flag1 = false;
    private int j;
    private double uOldSave;
    
    public MCMoveRotateMolecule3D(PotentialMaster potentialMaster, Space space) {
        super(potentialMaster, 1);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        rotationTensor = (Space3D.RotationTensor)space.makeRotationTensor();
        r0 = (Space3D.Vector)space.makeVector();
       
        setStepSizeMax(Math.PI);
        setStepSizeMin(0.0);
        setStepSize(Math.PI/2.0);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);
        //TODO enable this
        //set directive to exclude intramolecular contributions to the energy
//        iteratorDirective.addCriterion(new IteratorDirective.PotentialCriterion() {
//            public boolean excludes(Potential p) {return (p instanceof Potential1.Intramolecular);}
//        });
    }
     
    public boolean doTrial() {
        Phase phase = phases[0];
        if(phase.moleculeCount()==0) {molecule = null; return false;}
            
        molecule = phase.randomMolecule();
        energyMeter.setTarget(molecule);
        uOld = energyMeter.getDataAsScalar(phase);
        if(uOld < Double.MAX_VALUE) uOldSave = uOld;
        
        double dTheta = (2*Simulation.random.nextDouble() - 1.0)*stepSize;
        rotationTensor.setAxial(dTheta);
       
        leafAtomIterator.reset();
        r0.E(molecule.node.firstLeafAtom().coord.position());
        AtomActionTransform.doAction(leafAtomIterator, r0, rotationTensor);
            
        uNew = Double.NaN;
        return true;
    }//end of doTrial
    
    public double lnTrialRatio() {return 0.0;}
    
    public double lnProbabilityRatio() {
        energyMeter.setTarget(molecule);
        uNew = energyMeter.getDataAsScalar(phases[0]);
        if(uOld > Double.MAX_VALUE) uOld = uOldSave;
        return -(uNew - uOld)/temperature;
    }
    
    public void acceptNotify() {  /* do nothing */
            
    }
    
    public void rejectNotify() {
        leafAtomIterator.reset();
        
          rotationTensor.invert();
          leafAtomIterator.reset();
          AtomActionTransform.doAction(leafAtomIterator, r0, rotationTensor);
        
    }//end of rejectNotify
    
    public double energyChange(Phase phase) {
        if(this.phases[0] != phase) return 0.0;
        else return uNew - uOld;
    }

 
    public final AtomIterator affectedAtoms(Phase phase) {
        if(this.phases[0] != phase) return AtomIterator.NULL;
        affectedAtomIterator.setAtom(molecule);
        affectedAtomIterator.reset();
        return affectedAtomIterator;
    }


}
