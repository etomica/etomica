package etomica.association;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IAtomPositioned;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.ISimulation;
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
        double vOld = box.getBoundary().volume();
        uOld = energyMeter.getDataAsScalar();
        hOld = uOld + pressure*vOld;
        vScale = (2.*random.nextDouble()-1.)*stepSize;
        vNew = vOld * Math.exp(vScale); //Step in ln(V)
        double rScale = Math.exp(vScale/D);
        numAssociatedAtoms = associationManager.getAssociatedAtoms().getAtomCount();
//        System.out.println("doTrial, associated atoms = " +associationManager.getAssociatedAtoms());//list of the associated atoms
//        associationManager.initialize();
//        if (numAssociatedAtoms != associationManager.getAssociatedAtoms().getAtomCount()){//numAssociatedAtoms after move
//        	System.out.println("doTrial, associated atoms .....= " +associationManager.getAssociatedAtoms());//list of the associated atoms
//    		throw new RuntimeException("Oops");
//    	}
        if (numAssociatedAtoms %2 == 1){
        	IAtomList list = associationManager.getAssociatedAtoms();
        	for ( int i = 0; i< numAssociatedAtoms;i++){
        		System.out.println("list "+list.getAtom(i)+":"+ associationManager.getAssociatedAtoms(list.getAtom(i)));
        	}
//        	list = box.getLeafList();
//        	IAtomOriented a = (IAtomOriented)list.getAtom(354);
//        	System.out.println("atom 354 " +a.getPosition()+ " " +a.getOrientation().getDirection());
//        	a = (IAtomOriented)list.getAtom(68);
//        	System.out.println("atom 68 " +a.getPosition()+ " " +a.getOrientation().getDirection());
//        	a = (IAtomOriented)list.getAtom(409);
//        	System.out.println("atom 409 " +a.getPosition()+ " " +a.getOrientation().getDirection());
//        	System.out.println("bonded" +associationManager.getAssociationDefinition().isAssociated(list.getAtom(354), list.getAtom(68)));
//        	System.out.println("bonded" +associationManager.getAssociationDefinition().isAssociated(list.getAtom(354), list.getAtom(409)));
//        	System.out.println("bonded" +associationManager.getAssociationDefinition().isAssociated(list.getAtom(68), list.getAtom(409)));
//        	P2HardAssociationCone p = new P2HardAssociationCone(Space3D.getInstance(), 1.0, 1.0, 6.0, 16.0);
//        	p.setBox(box);
//        	System.out.println("energy 354-68 = "+p.energy(new AtomPair(list.getAtom(354), list.getAtom(68))));
//        	System.out.println("energy 354-409 = "+p.energy(new AtomPair(list.getAtom(354), list.getAtom(409))));
//        	System.out.println("energy 68-409 = "+p.energy(new AtomPair(list.getAtom(68), list.getAtom(409))));
        	throw new RuntimeException("***");
        }
        //System.out.println("rScale = "+rScale);
        scaleAtoms(rScale);//call the method
        uNew = energyMeter.getDataAsScalar();
        hNew = uNew + pressure*vNew;
        return true;
    }//end of doTrial
    
    protected void scaleAtoms(double rScale) {
    	IAtomList atomList = box.getLeafList();
    	for (int i=0; i< atomList.getAtomCount(); i++){
        	IAtom atom = atomList.getAtom(i);
        	IAtomList bondedAtoms = associationManager.getAssociatedAtoms(atom);
        	if (bondedAtoms.getAtomCount() == 0) {
        		continue;//go to another atom
        	}
        	if (atom.getLeafIndex() > bondedAtoms.getAtom(0).getLeafIndex()){
        		continue; //we skip the movement 
        	}
//        	if (atom.getLeafIndex() == 324 || atom.getLeafIndex() == 402 || bondedAtoms.getAtom(0).getLeafIndex() == 324 || bondedAtoms.getAtom(0).getLeafIndex() == 402) {
//    			System.out.println("atom1 = "+ atom+" "+((IAtomPositioned)atom).getPosition());
//    			System.out.println("bondedAtom = "+ bondedAtoms.getAtom(0)+" "+((IAtomPositioned)bondedAtoms.getAtom(0)).getPosition());
//    		}
        	r.Ev1Mv2(((IAtomPositioned)bondedAtoms.getAtom(0)).getPosition(),((IAtomPositioned)atom).getPosition());//position2 - position1
        	box.getBoundary().nearestImage(r);//choose the shorter distance
//        	if (atom.getLeafIndex() == 324 || atom.getLeafIndex() == 402 || bondedAtoms.getAtom(0).getLeafIndex() == 324 || bondedAtoms.getAtom(0).getLeafIndex() == 402) {
//        		System.out.println("rsquared "+ r.squared());
//        	}
        	((IAtomPositioned)bondedAtoms.getAtom(0)).getPosition().Ev1Pv2(((IAtomPositioned)atom).getPosition(), r);//move atom2 to the outside of box
        }
    	if (rScale > 1){//expanding
    		r.Ea1Tv1(rScale, box.getBoundary().getDimensions());
    		//System.out.println("box size = " +box.getBoundary().getDimensions());
    		//System.out.println("r = " +r);
            box.getBoundary().setDimensions(r);//scale the boundary
    	}
        
        for (int i=0; i< atomList.getAtomCount(); i++){
        	IAtom atom = atomList.getAtom(i);
        	IAtomList bondedAtoms = associationManager.getAssociatedAtoms(atom);
        	if (bondedAtoms.getAtomCount() == 0) {
//        		if (atom.getLeafIndex() == 324 || atom.getLeafIndex() == 402) {
//        			System.out.println("atom  = "+ atom +" " +((IAtomPositioned)atom).getPosition());
//        		}
        		((IAtomPositioned)atom).getPosition().TE(rScale);
//        		if (atom.getLeafIndex() == 324 || atom.getLeafIndex() == 402) {
//        			System.out.println("atom ... = "+ atom +" " +((IAtomPositioned)atom).getPosition());
//        		}
        		continue;//go to another atom
        	}
        	if (atom.getLeafIndex() > bondedAtoms.getAtom(0).getLeafIndex()){
        		continue; //we skip the movement 
        	}
//        	if (atom.getLeafIndex() == 324 || atom.getLeafIndex() == 402 || bondedAtoms.getAtom(0).getLeafIndex() == 324 || bondedAtoms.getAtom(0).getLeafIndex() == 402) {
//    			System.out.println("atom1 = "+ atom+" "+((IAtomPositioned)atom).getPosition());
//    			System.out.println("bondedAtom = "+ bondedAtoms.getAtom(0)+" "+((IAtomPositioned)bondedAtoms.getAtom(0)).getPosition());
//    		}
        	r.Ev1Mv2(((IAtomPositioned)bondedAtoms.getAtom(0)).getPosition(),((IAtomPositioned)atom).getPosition());//position2 - position1
        	box.getBoundary().nearestImage(r);//choose the shorter distance
//        	if (atom.getLeafIndex() == 324 || atom.getLeafIndex() == 402 || bondedAtoms.getAtom(0).getLeafIndex() == 324 || bondedAtoms.getAtom(0).getLeafIndex() == 402) {
//        		System.out.println("rsquared "+ r.squared());
//        	}
        	((IAtomPositioned)bondedAtoms.getAtom(0)).getPosition().Ev1Pv2(((IAtomPositioned)atom).getPosition(), r);//move atom2 to the outside of box
        	r.TE(0.5);//half of the separation distance
        	r.PE(((IAtomPositioned)atom).getPosition());//atom1 position + half of the separation distance
        	r.PE(box.getBoundary().centralImage(r));// position in the box, prevent the position outside of the box
        	r.TE(rScale-1);
        	
        	((IAtomPositioned)atom).getPosition().PE(r);//new position of atom1
        	((IAtomPositioned)bondedAtoms.getAtom(0)).getPosition().PE(r);//new position of atom2
        	//IVector dr = box.getBoundary().centralImage(((IAtomPositioned)bondedAtoms.getAtom(0)).getPosition());
        	//((IAtomPositioned)bondedAtoms.getAtom(0)).getPosition().PE(dr)
//        	if (atom.getLeafIndex() == 324 || atom.getLeafIndex() == 402 || bondedAtoms.getAtom(0).getLeafIndex() == 324 || bondedAtoms.getAtom(0).getLeafIndex() == 402) {
//    			System.out.println("atom1 .....= "+ atom+" "+((IAtomPositioned)atom).getPosition());
//    			System.out.println("bondedAtom .....= "+ bondedAtoms.getAtom(0)+" "+((IAtomPositioned)bondedAtoms.getAtom(0)).getPosition());
//    			r.Ev1Mv2(((IAtomPositioned)bondedAtoms.getAtom(0)).getPosition(),((IAtomPositioned)atom).getPosition());//position2 - position1
//            	box.getBoundary().nearestImage(r);//choose the shorter distance
//            	System.out.println("after rsquared "+ r.squared());
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
//    	int newnumAssociatedAtoms = associationManager.getAssociatedAtoms().getAtomCount();
//    	associationManager.initialize();
//    	if (newnumAssociatedAtoms != associationManager.getAssociatedAtoms().getAtomCount()){
//    		throw new RuntimeException("wrong");
//    	}
//    	System.out.println("getA, number of associated atoms = " +associationManager.getAssociatedAtoms().getAtomCount());
//    	System.out.println("number of num.associated atoms = " +numAssociatedAtoms);
//    	if (numAssociatedAtoms > associationManager.getAssociatedAtoms().getAtomCount()){
//    		System.out.println("getA, associated atoms = " +associationManager.getAssociatedAtoms());//list of the associated atoms
//    		throw new RuntimeException ("!!!");
//    	}
    	if (numAssociatedAtoms != associationManager.getAssociatedAtoms().getAtomCount()){//numAssociatedAtoms after move
    		return 0;
    	}
        return Math.exp((box.getMoleculeList().getMoleculeCount()-numAssociatedAtoms/2+1)*vScale);//num.monomer+num.dimer(not the total particle)
    }
    
    public double getB() {
        return -(hNew - hOld);
    }
    
    public void acceptNotify() {  /* do nothing */}
    
    public void rejectNotify() {
    	double rScale = Math.exp(vScale/D);
        scaleAtoms(1/rScale);
//        int newnumAssociatedAtoms = associationManager.getAssociatedAtoms().getAtomCount();
//    	associationManager.initialize();
//    	if (newnumAssociatedAtoms != associationManager.getAssociatedAtoms().getAtomCount()){
//    		System.out.println("newnumAssociatedAtoms= "+newnumAssociatedAtoms);
//    		System.out.println("associationManager.getAssociatedAtoms().getAtomCount()= "+associationManager.getAssociatedAtoms().getAtomCount());
//    		throw new RuntimeException("rejectNotify wrong");
//    	}
    }

    public double energyChange() {return uNew - uOld;}
    
    public AtomIterator affectedAtoms() {
        return affectedAtomIterator;
    }

    public void setPressure(double p) {pressure = p;}
    public final double getPressure() {return pressure;}
    public Dimension getPressureDimension() {return Pressure.DIMENSION;}
}