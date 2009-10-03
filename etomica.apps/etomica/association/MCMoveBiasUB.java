package etomica.association;
import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IRandom;
import etomica.api.IVectorMutable;
import etomica.atom.IAtomOriented;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.space.ISpace;
import etomica.space.IVectorRandom;

public class MCMoveBiasUB extends MCMoveBox {
    
    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private BiasVolume biasVolume;
    private IRandom random;
    private AssociationManager associationManager;
    private MeterPotentialEnergy meterPotentialEnergy;
    private int ni, Nai;
    private IAtom atomA;
    private IVectorRandom orientation;
    private boolean isbonding;
    private final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    private double uOld;
    private double uNew;
    private IVectorMutable oldPosition;
    private IVectorMutable oldDirection;

    public MCMoveBiasUB(PotentialMasterCell potentialMaster, BiasVolume bv, IRandom random, ISpace space) {
        super(potentialMaster);//variable
        biasVolume = bv;
        this.random =random;
        meterPotentialEnergy = new MeterPotentialEnergy(potentialMaster);
        orientation = (IVectorRandom)space.makeVector();
        oldPosition = space.makeVector();
        oldDirection = space.makeVector();
    }
    
    public void setBox(IBox box){
    	super.setBox(box);
    	meterPotentialEnergy.setBox(box);
    }
    
    public void setAssociationManager(AssociationManager associationManager){
    	this.associationManager = associationManager;
    }
    
    public boolean doTrial() {
        
        int N = box.getMoleculeList().getMoleculeCount();
        if(N < 2) return false;
        /*
        choose bonding or unbonding
        */
        isbonding = random.nextInt(2) == 0;
         
        if (isbonding) { // bonding
            IAtomList atoms = box.getLeafList();
            atomA = atoms.getAtom(random.nextInt(atoms.getAtomCount()));
            IAtom atomB = atomA;
            while(atomB == atomA) atomB = atoms.getAtom(random.nextInt(atoms.getAtomCount()));
            meterPotentialEnergy.setTarget(atomA);
            uOld = meterPotentialEnergy.getDataAsScalar();
            ni = associationManager.getAssociatedAtoms(atomA).getAtomCount();
            Nai = associationManager.getAssociatedAtoms().getAtomCount();
            oldPosition.E(atomA.getPosition());
            oldDirection.E(((IAtomOriented)atomA).getOrientation().getDirection());
            biasVolume.biasInsert(atomA, atomB);
            //System.out.println("atomA = " +atomA + " atomB = " +atomB);
         
        }//end bonding
        else { // unbonding,breaking bond
        	IAtomList atoms = associationManager.getAssociatedAtoms();
        	if (atoms.getAtomCount() == 0) {
        		return false;
        	}
            atomA = atoms.getAtom(random.nextInt(atoms.getAtomCount()));
            ni = associationManager.getAssociatedAtoms(atomA).getAtomCount();
            Nai = associationManager.getAssociatedAtoms().getAtomCount();
            meterPotentialEnergy.setTarget(atomA);
            uOld = meterPotentialEnergy.getDataAsScalar();
            oldPosition.E(atomA.getPosition());
            oldDirection.E(((IAtomOriented)atomA).getOrientation().getDirection());
            ((IVectorRandom)atomA.getPosition()).setRandomCube(random);
            atomA.getPosition().TE(box.getBoundary().getBoxSize());//translate
            orientation.setRandomSphere(random);
            ((IAtomOriented)atomA).getOrientation().setDirection(orientation);//orientation
        }//end unbonding
        affectedAtomIterator.setAtom(atomA);
        //System.out.println("MCMoveBiasUB = " +atomA);
        return true;
    }//end doTrial
    public double getB() {
    	uNew = meterPotentialEnergy.getDataAsScalar();
    	if (atomA.getParentGroup().getIndex() == 10 || atomA.getParentGroup().getIndex() == 452){
    		System.out.println("MCMoveBiasUB "+atomA);
        	System.out.println("uOld-uNew = "+(uOld-meterPotentialEnergy.getDataAsScalar()));
        	System.out.println("uOld = "+uOld);
        	System.out.println("isbonding = "+isbonding);
        }
    	return uOld - uNew;
    }
    public double getA() {
    	int Naj = associationManager.getAssociatedAtoms().getAtomCount();
    	int N = box.getMoleculeList().getMoleculeCount();
    	double phi = biasVolume.biasVolume()/box.getBoundary().volume();
    	//if (Naj == 0) System.out.println("A1 = " +(ni*Nai/((N-1)*phi)));
    	if(Naj == 0) return ni*Nai/((N-1)*phi);//acceptance criteria
    	int nj = associationManager.getAssociatedAtoms(atomA).getAtomCount();
    	//if (Nai == 0) System.out.println("A2 = "+((N-1)*phi/(nj*Naj)));
        if(Nai == 0) return (N-1)*phi/(nj*Naj);
        int deltaj = (nj == 0) ? 0 : 1;
        int deltai = (ni == 0) ? 0 : 1;
        //System.out.println("a =" +((N-1)*phi*deltaj/Naj + ni)/((N-1)*phi*deltai/Nai + nj));
        //System.out.println("ni = " +ni +" nj = " +nj);
        //System.out.println("isbonding = "+isbonding);
//        if (Double.isNaN(((N-1)*phi*deltaj/Naj + ni)/((N-1)*phi*deltai/Nai + nj))){
//        	throw new RuntimeException ("Oops");
//        }
        return ((N-1)*phi*deltaj/Naj + ni)/((N-1)*phi*deltai/Nai + nj);    
    }

	public AtomIterator affectedAtoms() {
		return affectedAtomIterator;
	}

	public double energyChange() {
		return uNew - uOld;
	}

	public void acceptNotify() {
		//System.out.println("accepted");
		
	}

	public void rejectNotify() {
		atomA.getPosition().E(oldPosition);
		((IAtomOriented)atomA).getOrientation().setDirection(oldDirection);
		//System.out.println("rejected");
	}
 
}//end of MCMoveBiasUB