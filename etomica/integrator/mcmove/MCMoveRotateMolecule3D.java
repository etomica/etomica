package etomica.integrator.mcmove;
import etomica.api.IVector;
import etomica.atom.AtomSet;
import etomica.atom.AtomSource;
import etomica.atom.AtomSourceRandomMolecule;
import etomica.atom.IAtomPositioned;
import etomica.atom.IMolecule;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.potential.PotentialMaster;
import etomica.space.RotationTensor;
import etomica.util.IRandom;


public class MCMoveRotateMolecule3D extends MCMoveBoxStep {
    
    private static final long serialVersionUID = 2L;
    private final MeterPotentialEnergy energyMeter;
    protected final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    protected final IRandom random;
    protected AtomSource moleculeSource;
    
    protected transient double uOld;
    protected transient double uNew = Double.NaN;
    protected transient IMolecule molecule;
    protected transient IVector r0;
    protected transient RotationTensor rotationTensor;
    public int count;
    public int count1;
    public boolean flag = false;
    public boolean flag1 = false;
    
    public MCMoveRotateMolecule3D(PotentialMaster potentialMaster, IRandom random) {
        super(potentialMaster);
        this.random = random;
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        rotationTensor = potentialMaster.getSpace().makeRotationTensor();
        r0 = potentialMaster.getSpace().makeVector();
       
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
//        System.out.println("doTrial MCMoveRotateMolecule called");
        
        if(box.moleculeCount()==0) {molecule = null; return false;}
            
        molecule = (IMolecule)moleculeSource.getAtom();
        energyMeter.setTarget(molecule);
        uOld = energyMeter.getDataAsScalar();
        if(Double.isInfinite(uOld)) {
            throw new RuntimeException("Overlap in initial state");
        }
        
        double dTheta = (2*random.nextDouble() - 1.0)*stepSize;
        rotationTensor.setAxial(r0.getD() == 3 ? random.nextInt(3) : 2,dTheta);

        r0.E(molecule.getType().getPositionDefinition().position(molecule));
        doTransform();
        
        energyMeter.setTarget(molecule);
        uNew = energyMeter.getDataAsScalar();
        return true;
    }//end of doTrial
    
    protected void doTransform() {
        AtomSet childList = molecule.getChildList();
        for (int iChild = 0; iChild<childList.getAtomCount(); iChild++) {
            IAtomPositioned a = (IAtomPositioned)childList.getAtom(iChild);
            IVector r = a.getPosition();
            r.ME(r0);
            box.getBoundary().nearestImage(r);
            rotationTensor.transform(r);
            r.PE(r0);
        }
    }
    
    public double getA() {return 1.0;}
    
    public double getB() {

        return -(uNew - uOld);
    }
    
    public void acceptNotify() {  /* do nothing */
            
    }
    
    public void rejectNotify() {
        rotationTensor.invert();
        doTransform();
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
     * @see etomica.integrator.MCMove#setBox(etomica.Box[])
     */
    public void setBox(Box p) {
        super.setBox(p);
        energyMeter.setBox(p);
        moleculeSource.setBox(p);
    }
}
