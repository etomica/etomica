/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;
import etomica.action.AtomActionTranslateBy;
import etomica.action.MoleculeChildAtomAction;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.util.random.IRandom;
import etomica.space.Vector;
import etomica.atom.MoleculeArrayList;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.atom.iterator.MoleculeIterator;
import etomica.atom.iterator.MoleculeIteratorSinglet;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.integrator.mcmove.MCMoveMolecular;
import etomica.models.OPLS.SpeciesAceticAcid;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space3d.RotationTensor3D;

public class MCMoveBiasUBMolecule extends MCMoveBox implements MCMoveMolecular{
    
    /** MC trials for molecules that implement the unbonding-bonding (UB) association bias algorithm 
     * for making and breaking an association with high efficiency
     * 
     * @author Hye Min Kim
	 */
	private static final long serialVersionUID = 1L;
	private BiasVolumeMolecule biasVolume;
    private IRandom random;
    private AssociationManagerMolecule associationManager;
    private IAssociationHelperMolecule associationHelper;
    private MeterPotentialEnergy meterPotentialEnergy;
    private int ni, Nai;
    private IMolecule moleculeA;
    private boolean isbonding;
    private final AtomIteratorArrayListSimple affectedAtomIterator = new AtomIteratorArrayListSimple();
    protected final MoleculeIteratorSinglet affectedMoleculeIterator = new MoleculeIteratorSinglet();
    private double uOld;
    private double uNew;
    private Vector[] oldPosition;
    protected final MoleculeArrayList smerList;
    protected final Space space;
    protected final RotationTensor3D rotationTensor;
    protected Vector groupTranslationVector;
    protected MoleculeChildAtomAction moveMoleculeAction;
    
    public MCMoveBiasUBMolecule(PotentialMaster potentialMaster, BiasVolumeMolecule bv, IRandom random, Space space) {
        super(potentialMaster);//variable
        biasVolume = bv;
        this.random =random;
        this.smerList = new MoleculeArrayList();
        this.space = space;
        meterPotentialEnergy = new MeterPotentialEnergy(potentialMaster);
        rotationTensor = (RotationTensor3D)(space.makeRotationTensor());
        AtomActionTranslateBy translator = new AtomActionTranslateBy(space);
        groupTranslationVector = translator.getTranslationVector();
        moveMoleculeAction = new MoleculeChildAtomAction(translator);
        perParticleFrequency = true;// the frequency of the move is increasing with the system size
        oldPosition = new Vector[0];
    }
    
    public void setBox(Box box){
    	super.setBox(box);
    	meterPotentialEnergy.setBox(box);
    }
    
    public void setAssociationManager(AssociationManagerMolecule associationManager, IAssociationHelperMolecule associationHelper){
    	this.associationManager = associationManager;
    	this.associationHelper = associationHelper;
    }
    
    public boolean doTrial() {//asking the biasVolume class form a bond or break bond
        int N = box.getMoleculeList().getMoleculeCount();
        if(N < 2) return false;
        /*
        choose bonding or unbonding
        */
        isbonding = random.nextInt(2) == 0;
         
        if (isbonding) { // bonding
        	
            IMoleculeList molecules = box.getMoleculeList();
            moleculeA = molecules.getMolecule(random.nextInt(molecules.getMoleculeCount()));
            IAtomList atoms = moleculeA.getChildList();
            if (oldPosition.length < atoms.getAtomCount()){
            	oldPosition = new Vector[atoms.getAtomCount()];
                for (int i = 0; i<atoms.getAtomCount();i+=1){
                	oldPosition[i] = space.makeVector();
                }
            }
            IMolecule moleculeB = moleculeA;
            while(moleculeB == moleculeA) moleculeB = molecules.getMolecule(random.nextInt(molecules.getMoleculeCount()));
            meterPotentialEnergy.setTarget(moleculeA);
            uOld = meterPotentialEnergy.getDataAsScalar();
            ni = associationManager.getAssociatedMolecules(moleculeA).getMoleculeCount();
            Nai = associationManager.getAssociatedMolecules().getMoleculeCount();
            for (int i = 0; i<atoms.getAtomCount();i+=1){
            	oldPosition[i].E(atoms.getAtom(i).getPosition());
            }
            biasVolume.biasInsert(moleculeA, moleculeB);
        }//end bonding
        else { // unbonding,breaking bond
        	IMoleculeList molecules = associationManager.getAssociatedMolecules();//associated molecules
        	if (molecules.getMoleculeCount() == 0) {
        		return false;
        	}
            moleculeA = molecules.getMolecule(random.nextInt(molecules.getMoleculeCount()));
            IAtomList atoms = moleculeA.getChildList();
            if (oldPosition.length < atoms.getAtomCount()){
            	oldPosition = new Vector[atoms.getAtomCount()];
                for (int i = 0; i<atoms.getAtomCount();i+=1){
                	oldPosition[i] = space.makeVector();
                }
            }
            ni = associationManager.getAssociatedMolecules(moleculeA).getMoleculeCount();
            Nai = associationManager.getAssociatedMolecules().getMoleculeCount();
            meterPotentialEnergy.setTarget(moleculeA);
            uOld = meterPotentialEnergy.getDataAsScalar();
            for (int i = 0; i<atoms.getAtomCount();i+=1){
            	oldPosition[i].E(atoms.getAtom(i).getPosition());
            }
            doTransform(moleculeA, oldPosition[SpeciesAceticAcid.indexC]);//random rotation of moleculeA
            groupTranslationVector.setRandomCube(random);
            groupTranslationVector.TE(box.getBoundary().getBoxSize());//translate
            moveMoleculeAction.actionPerformed(moleculeA);
        }//end unbonding
        affectedAtomIterator.setList(moleculeA.getChildList());
        return true;
    }//end doTrial
    public double getB() {
    	uNew = meterPotentialEnergy.getDataAsScalar();
    	return uOld - uNew;
    }
    public double getA() {
    	int Naj = associationManager.getAssociatedMolecules().getMoleculeCount();
    	int N = box.getMoleculeList().getMoleculeCount();
    	double phi = biasVolume.biasVolume()/box.getBoundary().volume()*N;

    	if(Naj == 0) return ni*Nai/((N-1)*phi);//acceptance criteria
    	int nj = associationManager.getAssociatedMolecules(moleculeA).getMoleculeCount();
        if(Nai == 0) return (N-1)*phi/(nj*Naj);
        int deltaj = (nj == 0) ? 0 : 1;
        int deltai = (ni == 0) ? 0 : 1;
        if (associationHelper.populateList(smerList,moleculeA,true)){
        	return 0;
        }
        return ((N-1)*phi*deltaj/Naj + ni)/((N-1)*phi*deltai/Nai + nj);    
    }
    
   	public AtomIterator affectedAtoms() {
		return affectedAtomIterator;
	}
   	
    protected void doTransform(IMolecule molecule, Vector r0) {
        IAtomList childList = molecule.getChildList();
        rotationTensor.setAxial(random.nextInt(3), random.nextDouble()*2*Math.PI);//axis for rotation is x, y or z, and the amount of rotation is up to 2pi
        for (int iChild = 0; iChild<childList.getAtomCount(); iChild++) {
            IAtom a = childList.getAtom(iChild);
            Vector r = a.getPosition();
            r.ME(r0);
            rotationTensor.transform(r);
            r.PE(r0);
        }
    }

	public double energyChange() {
		return uNew - uOld;
	}

	public void acceptNotify() {
	}

	public void rejectNotify() {
        for (int i = 0; i<moleculeA.getChildList().getAtomCount();i+=1){
        	moleculeA.getChildList().getAtom(i).getPosition().E(oldPosition[i]);
        }
	}

	public MoleculeIterator affectedMolecules(Box box) {
        affectedMoleculeIterator.setMolecule(moleculeA);
        return affectedMoleculeIterator;
	}
}
