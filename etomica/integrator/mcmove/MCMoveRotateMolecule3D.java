package etomica.integrator.mcmove;
import etomica.action.AtomTransform;
import etomica.atom.Atom;
import etomica.atom.AtomSource;
import etomica.atom.AtomSourceRandomMolecule;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.atom.iterator.AtomIteratorTree;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Vector3D;
import etomica.util.IRandom;


public class MCMoveRotateMolecule3D extends MCMovePhaseStep {
    
    private static final long serialVersionUID = 2L;
    private final MeterPotentialEnergy energyMeter;
    protected final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    private final IRandom random;
    protected AtomSource moleculeSource;
    
    protected final AtomIteratorTree leafAtomIterator = new AtomIteratorTree();
    
    protected transient double uOld;
    protected transient double uNew = Double.NaN;
    protected transient Atom molecule;
    protected transient Vector3D r0;
    protected transient RotationTensor3D rotationTensor;
    public int count;
    public int count1;
    public boolean flag = false;
    public boolean flag1 = false;
    private double uOldSave;
    
    public MCMoveRotateMolecule3D(Simulation sim) {
        this(sim.getPotentialMaster(), sim.getRandom());
    }
    
    public MCMoveRotateMolecule3D(PotentialMaster potentialMaster, IRandom random) {
        super(potentialMaster);
        this.random = random;
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        rotationTensor = (RotationTensor3D)potentialMaster.getSpace().makeRotationTensor();
        r0 = (Vector3D)potentialMaster.getSpace().makeVector();
       
        setStepSizeMax(Math.PI);
        setStepSizeMin(0.0);
        setStepSize(Math.PI/2.0);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);
        moleculeSource = new AtomSourceRandomMolecule();
        ((AtomSourceRandomMolecule)moleculeSource).setRandom(random);
        //TODO enable this
        //set directive to exclude intramolecular contributions to the energy
//        iteratorDirective.addCriterion(new IteratorDirective.PotentialCriterion() {
//            public boolean excludes(Potential p) {return (p instanceof Potential1.Intramolecular);}
//        });
    }
     
    public boolean doTrial() {
        if(phase.moleculeCount()==0) {molecule = null; return false;}
            
        molecule = moleculeSource.getAtom();
        energyMeter.setTarget(molecule);
        uOld = energyMeter.getDataAsScalar();
        if(uOld < Double.MAX_VALUE) uOldSave = uOld;
        
        double dTheta = (2*random.nextDouble() - 1.0)*stepSize;
        rotationTensor.setAxial(random.nextInt(3),dTheta);

        leafAtomIterator.setRoot(molecule);
        leafAtomIterator.reset();
        r0.E(molecule.getType().getPositionDefinition().position(molecule));
        AtomTransform.doTransform(leafAtomIterator, r0, rotationTensor);
            
        uNew = Double.NaN;
        return true;
    }//end of doTrial
    
    public double getA() {return 1.0;}
    
    public double getB() {
        energyMeter.setTarget(molecule);
        uNew = energyMeter.getDataAsScalar();
        if(uOld > Double.MAX_VALUE) uOld = uOldSave;
        return -(uNew - uOld);
    }
    
    public void acceptNotify() {  /* do nothing */
            
    }
    
    public void rejectNotify() {
        rotationTensor.invert();
        leafAtomIterator.reset();
        AtomTransform.doTransform(leafAtomIterator, r0, rotationTensor);
    }
    
    public double energyChange() {
        return uNew - uOld;
    }

 
    public final AtomIterator affectedAtoms() {
        affectedAtomIterator.setAtom(molecule);
        affectedAtomIterator.reset();
        return affectedAtomIterator;
    }


    /* (non-Javadoc)
     * @see etomica.integrator.MCMove#setPhase(etomica.Phase[])
     */
    public void setPhase(Phase p) {
        super.setPhase(p);
        energyMeter.setPhase(p);
    }
}
