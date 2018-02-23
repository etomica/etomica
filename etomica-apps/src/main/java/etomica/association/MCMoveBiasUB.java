/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

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
    private int maxLength = Integer.MAX_VALUE;
    private IAtom atomA;
    private Vector orientation;
    private boolean isbonding;
    private final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    private double uOld;
    private double uNew;
    private Vector oldPosition;
    private Vector oldDirection;
    protected final AtomArrayList smerList;
    protected final Vector dr;

    public MCMoveBiasUB(PotentialMasterCell potentialMaster, BiasVolume bv, IRandom random, Space space) {
        super(potentialMaster);//variable
        biasVolume = bv;
        this.random =random;
        this.smerList = new AtomArrayList();
        this.dr = space.makeVector();
        meterPotentialEnergy = new MeterPotentialEnergy(potentialMaster);
        orientation = space.makeVector();
        oldPosition = space.makeVector();
        oldDirection = space.makeVector();
        perParticleFrequency = true;// the frequency of the move is increasing with the system size
    }
    
    public void setBox(Box box){
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
            atomA = atoms.get(random.nextInt(atoms.size()));
            IAtom atomB = atomA;
            while(atomB == atomA) atomB = atoms.get(random.nextInt(atoms.size()));
            meterPotentialEnergy.setTarget(atomA);
            uOld = meterPotentialEnergy.getDataAsScalar();
            ni = associationManager.getAssociatedAtoms(atomA).size();
            Nai = associationManager.getAssociatedAtoms().size();
            oldPosition.E(atomA.getPosition());
            oldDirection.E(((IAtomOriented)atomA).getOrientation().getDirection());
            biasVolume.biasInsert(atomA, atomB);
            //System.out.println("atomA = " +atomA + " atomB = " +atomB);
         
        }//end bonding
        else { // unbonding,breaking bond
        	IAtomList atoms = associationManager.getAssociatedAtoms();//associated atoms
        	if (atoms.size() == 0) {
        		return false;
        	}
            atomA = atoms.get(random.nextInt(atoms.size()));
            ni = associationManager.getAssociatedAtoms(atomA).size();
            Nai = associationManager.getAssociatedAtoms().size();
            meterPotentialEnergy.setTarget(atomA);
            uOld = meterPotentialEnergy.getDataAsScalar();
            oldPosition.E(atomA.getPosition());
            oldDirection.E(((IAtomOriented)atomA).getOrientation().getDirection());
            atomA.getPosition().setRandomCube(random);
            atomA.getPosition().TE(box.getBoundary().getBoxSize());//translate
            orientation.setRandomSphere(random);
            ((IAtomOriented)atomA).getOrientation().setDirection(orientation);//orientation
        }//end unbonding
        affectedAtomIterator.setAtom(atomA);
        //System.out.println("MCMoveBiasUB = " +atomA);
        return true;
    }//end doTrial
    public double getB() {
//    	if (atomA.getParentGroup().getIndex() == 10 || atomA.getParentGroup().getIndex() == 452){
//    		System.out.println("MCMoveBiasUB "+atomA);
//        	System.out.println("uOld-uNew = "+(uOld-meterPotentialEnergy.getDataAsScalar()));
//        	System.out.println("uOld = "+uOld);
//        	System.out.println("isbonding = "+isbonding);
//        }
    	return uOld - uNew;
    }

    public double getChi(double temperature) {
//    	System.out.print("isbonding= "+isbonding +" ");
//    	if (isbonding){
//    		System.out.print("\t\t" );
//    	}
//        if (associationManager.getAssociatedAtoms(atomA).getAtomCount() > 1) {
//        	return 0;
//        } 
//        if (associationManager.getAssociatedAtoms(atomA).getAtomCount() == 1){
//        	IAtom atomj = associationManager.getAssociatedAtoms(atomA).getAtom(0);
//        	if(associationManager.getAssociatedAtoms(atomj).getAtomCount() > 1){
//        		return 0;
//        	} 
//        }
//        else if(isbonding){
//        	throw new RuntimeException("wrong!!!");
//        }
    	int Naj = associationManager.getAssociatedAtoms().size();
    	int N = box.getMoleculeList().getMoleculeCount();
    	double phi = biasVolume.biasVolume()/box.getBoundary().volume()*N;
    	//if (Naj == 0) System.out.println("A1 = " +(ni*Nai/((N-1)*phi)));
    	if(Naj == 0) return ni*Nai/((N-1)*phi);//acceptance criteria
    	int nj = associationManager.getAssociatedAtoms(atomA).size();
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
        if (populateList(smerList) == 0){
        	return 0;
        }
        if (smerList.size() > maxLength) {
    		return 0.0;
		}
        uNew = meterPotentialEnergy.getDataAsScalar();
        return ((N - 1) * phi * deltaj / Naj + ni) / ((N - 1) * phi * deltai / Nai + nj) * Math.exp(-(uNew - uOld) / temperature);
    }
    
    public void setMaxLength(int i){
    	maxLength = i;
    }
    
    protected int populateList(AtomArrayList mySmerList){
    	mySmerList.clear();
    	mySmerList.add(atomA);
    	IAtomList bondList = associationManager.getAssociatedAtoms(atomA);
    	if (bondList.size() > 2){
    		return 0;
    	}
    	if (bondList.size() == 2){
    		IAtom atom0 = bondList.get(0);
    		IAtom atom1 = bondList.get(1);
    		dr.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());//dr = distance from the atom0 to atom1
        	box.getBoundary().nearestImage(dr);
        	double innerRadius = 0.8;
        	double minDistance = 2*(innerRadius*innerRadius)*(1+Math.cos(etomica.units.Degree.UNIT.toSim(27.0)));
        	if (dr.squared() < minDistance){
        		return 0;
        	}
    	}
    	if (bondList.size() == 0){
    		return 1;
    	}
    	IAtom thisAtom = bondList.get(0);
    	mySmerList.add(thisAtom);
    	IAtomList bondList1 = associationManager.getAssociatedAtoms(thisAtom);
    	if (bondList1.size() > 2){
    		return 0;
    	}
    	if (bondList1.size() == 2){
    		IAtom atom0 = bondList1.get(0);
    		IAtom atom1 = bondList1.get(1);
    		dr.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());//dr = distance from the atom0 to atom1
        	box.getBoundary().nearestImage(dr);
        	double innerRadius = 0.8;
        	double minDistance = 2*(innerRadius*innerRadius)*(1+Math.cos(etomica.units.Degree.UNIT.toSim(27.0)));
        	if (dr.squared() < minDistance){
        		return 0;
        	}
    	}
    	IAtom previousAtom = atomA;
    	while (bondList1.size() > 1){
    		IAtom nextAtom = bondList1.get(0);
    		if (nextAtom == previousAtom){
    			nextAtom = bondList1.get(1);
    		} 
    		if (nextAtom == atomA){
    			return 1;
    		}
    		mySmerList.add(nextAtom);
    		bondList1 = associationManager.getAssociatedAtoms(nextAtom);
    		if (bondList1.size() > 2){
        		return 0;
        	}
    		if (bondList1.size() == 2){
        		IAtom atom0 = bondList1.get(0);
        		IAtom atom1 = bondList1.get(1);
        		dr.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());//dr = distance from the atom0 to atom1
            	box.getBoundary().nearestImage(dr);
            	double innerRadius = 0.8;
            	double minDistance = 2*(innerRadius*innerRadius)*(1+Math.cos(etomica.units.Degree.UNIT.toSim(27.0)));
            	if (dr.squared() < minDistance){
            		return 0;
            	}
        	}
    		previousAtom = thisAtom;
    		thisAtom = nextAtom;
    	}
    	if (bondList.size()>1){
    		thisAtom = bondList.get(1);
        	mySmerList.add(thisAtom);
        	bondList1 = associationManager.getAssociatedAtoms(thisAtom);
        	if (bondList1.size() > 2){
        		return 0;
        	}
        	if (bondList1.size() == 2){
        		IAtom atom0 = bondList1.get(0);
        		IAtom atom1 = bondList1.get(1);
        		dr.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());//dr = distance from the atom0 to atom1
            	box.getBoundary().nearestImage(dr);
            	double innerRadius = 0.8;
            	double minDistance = 2*(innerRadius*innerRadius)*(1+Math.cos(etomica.units.Degree.UNIT.toSim(27.0)));
            	if (dr.squared() < minDistance){
            		return 0;
            	}
        	}
        	previousAtom = atomA;
        	while (bondList1.size() > 1){
        		IAtom nextAtom = bondList1.get(0);
        		if (nextAtom == previousAtom){
        			nextAtom = bondList1.get(1);
        		} 
        		mySmerList.add(nextAtom);
        		bondList1 = associationManager.getAssociatedAtoms(nextAtom);
        		if (bondList1.size() > 2){
            		return 0;
            	}
        		if (bondList1.size() == 2){
            		IAtom atom0 = bondList1.get(0);
            		IAtom atom1 = bondList1.get(1);
            		dr.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());//dr = distance from the atom0 to atom1
                	box.getBoundary().nearestImage(dr);
                	double innerRadius = 0.8;
                	double minDistance = 2*(innerRadius*innerRadius)*(1+Math.cos(etomica.units.Degree.UNIT.toSim(27.0)));
                	if (dr.squared() < minDistance){
                		return 0;
                	}
            	}
        		previousAtom = thisAtom;
        		thisAtom = nextAtom;
        	}
    	}
    	return 1;
    }

	public AtomIterator affectedAtoms() {
		return affectedAtomIterator;
	}

	public double energyChange() {
		return uNew - uOld;
	}

	public void acceptNotify() {
//		if (atomA.getLeafIndex()== 388 ||atomA.getLeafIndex()== 115 ){
//        	System.out.println("accepted UB moving atomA "+atomA);
//        	System.out.println("position 388 "+((IAtomPositioned)box.getLeafList().getAtom(388)).getPosition()+"position 115 "+((IAtomPositioned)box.getLeafList().getAtom(115)).getPosition());
//        }
		
	}

	public void rejectNotify() {
		atomA.getPosition().E(oldPosition);
		((IAtomOriented)atomA).getOrientation().setDirection(oldDirection);
//		if (atomA.getLeafIndex()== 388 ||atomA.getLeafIndex()== 115 ){
//        	System.out.println("rejected UB moving atomA "+atomA);
//        	System.out.println("position 388 "+((IAtomPositioned)box.getLeafList().getAtom(388)).getPosition()+"position 115 "+((IAtomPositioned)box.getLeafList().getAtom(115)).getPosition());
//        }
	}
}
