package etomica.association;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IAtomPositioned;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.space.ISpace;
import etomica.units.Dimension;
import etomica.units.Pressure;

/**
 * Standard Monte Carlo volume-change move for simulations in the NPT ensemble.
 *
 * @author David Kofke
 */
public class MCMoveVolumeAssociated extends MCMoveBoxStep {
    
    private static final long serialVersionUID = 2L;
    protected double pressure;
    private MeterPotentialEnergy energyMeter;
    private final int D;
    private IRandom random;
    protected final AtomIteratorLeafAtoms affectedAtomIterator;
    protected final IVectorMutable r;
    private transient double uOld, hOld, vNew, vScale, hNew;
    private transient double uNew = Double.NaN;
    protected AssociationManager associationManager;
    protected int numAssociatedAtoms;

    public MCMoveVolumeAssociated(ISimulation sim, IPotentialMaster potentialMaster,
    		            ISpace _space) {
        this(potentialMaster, sim.getRandom(), _space, 1.0);
    }
    
    /**
     * @param potentialMaster an appropriate PotentialMaster instance for calculating energies
     * @param space the governing space for the simulation
     */
    public MCMoveVolumeAssociated(IPotentialMaster potentialMaster, IRandom random,
    		            ISpace _space, double pressure) {
        super(potentialMaster);
        this.random = random;
        this.D = _space.D();
        r = _space.makeVector();
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.10);
        setPressure(pressure);
        energyMeter.setIncludeLrc(true);
        affectedAtomIterator = new AtomIteratorLeafAtoms();
    }
    
    public void setBox(IBox p) {
        super.setBox(p);
        energyMeter.setBox(p);
        affectedAtomIterator.setBox(p);
    }
    
    public boolean doTrial() {
    	//System.out.println("doTrial, number of associated atoms = " +associationManager.getAssociatedAtoms().getAtomCount());
        double vOld = box.getBoundary().volume();
        uOld = energyMeter.getDataAsScalar();
        hOld = uOld + pressure*vOld;
        vScale = (2.*random.nextDouble()-1.)*stepSize;
        vNew = vOld * Math.exp(vScale); //Step in ln(V)
        double rScale = Math.exp(vScale/D);
        numAssociatedAtoms = associationManager.getAssociatedAtoms().getAtomCount();
        scaleAtoms(rScale);//call the method
        uNew = energyMeter.getDataAsScalar();
        hNew = uNew + pressure*vNew;
        return true;
    }//end of doTrial
    
    protected void scaleAtoms(double rScale) {
    	if (rScale > 1){//expanding
    		r.Ea1Tv1(rScale, box.getBoundary().getDimensions());
    		//System.out.println("box size = " +box.getBoundary().getDimensions());
    		//System.out.println("r = " +r);
            box.getBoundary().setDimensions(r);//scale the boundary
    	}
        IAtomList atomList = box.getLeafList();
        for (int i=0; i< atomList.getAtomCount(); i++){
        	IAtom atom = atomList.getAtom(i);
        	IAtomList bondedAtoms = associationManager.getAssociatedAtoms(atom);
        	if (bondedAtoms.getAtomCount() == 0) {
//        		if (atom.getLeafIndex() == 501) {
//        			System.out.println("atom = "+ ((IAtomPositioned)atom).getPosition());
//        		}
        		((IAtomPositioned)atom).getPosition().TE(rScale);
//        		if (atom.getLeafIndex() == 501) {
//        			System.out.println("atom.. = "+ ((IAtomPositioned)atom).getPosition());
//        		}
        		continue;//go to another atom
        	}
        	if (atom.getLeafIndex() > bondedAtoms.getAtom(0).getLeafIndex()){
        		continue; //we skip the movement 
        	}
//        	if (atom.getLeafIndex() == 501) {
//    			System.out.println("atom1 = "+ ((IAtomPositioned)atom).getPosition());
//    		}
//        	if (bondedAtoms.getAtom(0).getLeafIndex() == 501) {
//    			System.out.println("bondedAtom = "+ ((IAtomPositioned)bondedAtoms.getAtom(0)).getPosition());
//    		}
        	r.Ev1Mv2(((IAtomPositioned)bondedAtoms.getAtom(0)).getPosition(),((IAtomPositioned)atom).getPosition());//position2 - position1
        	box.getBoundary().nearestImage(r);//choose the shorter distance
        	((IAtomPositioned)bondedAtoms.getAtom(0)).getPosition().Ev1Pv2(((IAtomPositioned)atom).getPosition(), r);//move atom2 to the outside of box
        	r.TE(0.5);//half of the separation distance
        	r.PE(((IAtomPositioned)atom).getPosition());//atom1 position + half of the separation distance
        	r.PE(box.getBoundary().centralImage(r));// position in the box, prevent the position outside of the box
        	r.TE(1-rScale);
        	
        	((IAtomPositioned)atom).getPosition().PE(r);//new position of atom1
        	((IAtomPositioned)bondedAtoms.getAtom(0)).getPosition().PE(r);//new position of atom2
        	IVector dr = box.getBoundary().centralImage(((IAtomPositioned)bondedAtoms.getAtom(0)).getPosition());
        	((IAtomPositioned)bondedAtoms.getAtom(0)).getPosition().PE(dr);
//        	if (atom.getLeafIndex() == 501) {
//    			System.out.println("atom1... = "+ ((IAtomPositioned)atom).getPosition());
//    		}
//        	if (bondedAtoms.getAtom(0).getLeafIndex() == 501) {
//    			System.out.println("bondedAtom... = "+ ((IAtomPositioned)bondedAtoms.getAtom(0)).getPosition());
//    		}
        }
        if (rScale < 1){//compressing
    		r.Ea1Tv1(rScale, box.getBoundary().getDimensions());
    		//System.out.println("box size = " +box.getBoundary().getDimensions());
    		//System.out.println("r = " +r);
            box.getBoundary().setDimensions(r);//scale the boundary
    	}
    }
    public void setAssociationManager(AssociationManager associationManager) {
    	this.associationManager = associationManager;
    }
    
    public double getA() {
    	System.out.println("getA, number of associated atoms = " +associationManager.getAssociatedAtoms().getAtomCount());
    	if (numAssociatedAtoms != associationManager.getAssociatedAtoms().getAtomCount()){//numAssociatedAtoms after move
    		return 0;
    	}
        return Math.exp((box.getMoleculeList().getMoleculeCount()+1)*vScale);
    }
    
    public double getB() {
        return -(hNew - hOld);
    }
    
    public void acceptNotify() {  /* do nothing */}
    
    public void rejectNotify() {
    	double rScale = Math.exp(vScale/D);
        scaleAtoms(1/rScale);
    }

    public double energyChange() {return uNew - uOld;}
    
    public AtomIterator affectedAtoms() {
        return affectedAtomIterator;
    }

    public void setPressure(double p) {pressure = p;}
    public final double getPressure() {return pressure;}
    public Dimension getPressureDimension() {return Pressure.DIMENSION;}
}