package etomica;
import etomica.action.AtomActionTransform;


public class MCMoveRotateMolecule3D extends MCMove {
    
    private final IteratorDirective iteratorDirective = new IteratorDirective(IteratorDirective.BOTH);
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
    
    public MCMoveRotateMolecule3D(IntegratorMC parentIntegrator) {
        super(parentIntegrator);
        rotationTensor = (Space3D.RotationTensor)parentIntegrator.space.makeRotationTensor();
        r0 = (Space3D.Vector)parentIntegrator.space.makeVector();
       
        setStepSizeMax(Math.PI);
        setStepSizeMin(0.0);
        setStepSize(Math.PI/2.0);
        setPerParticleFrequency(true);
        iteratorDirective.includeLrc = false;
        //set directive to exclude intramolecular contributions to the energy
        iteratorDirective.addCriterion(new IteratorDirective.PotentialCriterion() {
            public boolean excludes(Potential p) {return (p instanceof Potential1.Intramolecular);}
        });
    }
     
    public boolean doTrial() {
        if(phase.moleculeCount()==0) {molecule = null; return false;}
            
        molecule = phase.randomMolecule();
        leafAtomIterator.setRoot(molecule);
		uOld = potential.calculate(phase, iteratorDirective.set(molecule), energy.reset()).sum();
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
		uNew = potential.calculate(phase, iteratorDirective.set(molecule), energy.reset()).sum();//not thread safe for multiphase systems
        if(uOld > Double.MAX_VALUE) uOld = uOldSave;
        return -(uNew - uOld)/parentIntegrator.temperature;
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
        if(this.phase != phase) return 0.0;
        else return uNew - uOld;
    }

 
    public final AtomIterator affectedAtoms(Phase phase) {
        if(this.phase != phase) return AtomIterator.NULL;
        affectedAtomIterator.setBasis(molecule);
        affectedAtomIterator.reset();
        return affectedAtomIterator;
    }


}
