package etomica.models.nitrogen;

import etomica.action.AtomActionTranslateBy;
import etomica.action.MoleculeChildAtomAction;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IPotentialMaster;
import etomica.api.IPotentialMolecular;
import etomica.api.IRandom;
import etomica.atom.AtomArrayList;
import etomica.atom.MoleculePair;
import etomica.atom.MoleculeSource;
import etomica.atom.MoleculeSourceRandomMolecule;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.space.ISpace;
import etomica.space.IVectorRandom;

/**
 * Standard Monte Carlo molecule-displacement trial move.  Two molecules are moved at a
 * time in such a way that the geometric center of the system is not changed.
 *
 * @author Tai Boon Tan
 */
public class MCMoveMoleculeCoupledSuperBox extends MCMoveBoxStep {

    
    public MCMoveMoleculeCoupledSuperBox(IPotentialMaster potentialMaster, IRandom nRandom,
    		                     ISpace _space, IBox box, int nC, int basis, CoordinateDefinitionNitrogenSuperBox coordinateDef){
        super(potentialMaster);
        this.random = nRandom;
        this.box = box;
        this.basisCell = coordinateDef.getBasisCells();
        
        super.setBox(box);
        
        moleculeSource = new MoleculeSourceRandomMolecule();
        ((MoleculeSourceRandomMolecule)moleculeSource).setRandom(random);
        moleculeSource.setBox(box);
        
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        energyMeter.setBox(box);
        
        affectedMoleculeList = new AtomArrayList();
        affectedMoleculeIterator = new AtomIteratorArrayListSimple(affectedMoleculeList);
        
        singleAction = new AtomActionTranslateBy(_space);
        groupTransVect = (IVectorRandom)singleAction.getTranslationVector();
        
        moveMoleculeAction = new MoleculeChildAtomAction(singleAction);
        
        pair = new MoleculePair();
        pairAB = new MoleculePair();
        
        perParticleFrequency = true;
        
        nA = (nC*nC*nC)*basis;
        molIndex = new int[27][nA];
        
        int layerNum = (nC*basis*3)*(3*nC);
        int oneThirdLayerNum = (nC*basis*3)*(nC);
        
        int axisNum = nC*basis*3;
        int cellAxisNum = nC*basis;
        int ix= 0;
        int iy= 0;
        int iz= 0;
        
        for (int iCell=0; iCell<molIndex.length; iCell++){
        	
        	int counter = 0;
        	int cellConst = 0;
        	
        	if(iCell%27 > 17){
        		ix = 2;
        	} else if (iCell%27 > 8 && iCell%27 <=17){
        		ix = 1;
        	} else {
        		ix = 0;
        	}
        	
        	if(iCell%9 > 5){
        		iy = 2;
        	} else if (iCell%9 > 2 && iCell%9 <=5){
        		iy = 1;
        	} else {
        		iy = 0;
        	}
        	
        	iz = iCell%3;
        
        	cellConst = iz*(cellAxisNum) + iy*oneThirdLayerNum + ix*nC*layerNum;
        	
	        for(int xnC=0; xnC<nC; xnC++){
		        for(int ynC=0; ynC<nC; ynC++){
			        for(int i=0; i<nC*basis; i++){
			        	molIndex[iCell][counter] = (i+(ynC*axisNum)+(xnC*layerNum)+cellConst);
			        	++counter;
			        }
		        }
	        }
        }
        
    }

    public void setBox(IBox newBox) {
 
    }
    
    public void setPotential(IPotentialMolecular newPotential){
        potential = newPotential;
    }
    
    public AtomIterator affectedAtoms() {
        affectedMoleculeList.clear();
        affectedMoleculeList.addAll(molecule0.getChildList());
        affectedMoleculeList.addAll(molecule1.getChildList());
        return affectedMoleculeIterator;
    }

    public double energyChange() {return uNew - uOld;}

    public void acceptNotify() {
        // I do believe nothing needs to happen here.
    }

    public boolean doTrial() {
//        System.out.println("doTrial MCMoveMoleculeCoupled called");
        
        randomMol0 = random.nextInt(nA);
        randomMol1 = random.nextInt(nA);
        //System.out.println("randomMol: " + randomMol0 + " " + randomMol1);
        molecule0 = basisCell[0].molecules.getMolecule(molIndex[13][randomMol0]);
        molecule1 = basisCell[0].molecules.getMolecule(molIndex[13][randomMol1]);
        
        if(molecule0==null || molecule1==null || molecule0==molecule1) return false;
        
        energyMeter.setTarget(molecule0);
        uOld = energyMeter.getDataAsScalar();
        
        //molSpeciesA
        
        double uCorrect = 0;
        for(int i=0; i<molIndex[0].length; i++){
        	molSpeciesA = basisCell[0].molecules.getMolecule(molIndex[13][i]);
        	pairAB.atom0 = molSpeciesA;
        	
        	if(molSpeciesA == molecule0 || molSpeciesA == molecule1) continue;
        	
        	for(int nCellBox=0; nCellBox<molIndex.length; nCellBox++){
        		if(nCellBox != 13){
        			molSpeciesB = basisCell[0].molecules.getMolecule(molIndex[nCellBox][randomMol0]);
        			pairAB.atom1 = molSpeciesB;
        			uCorrect += potential.energy(pairAB)/2;
        			
        			molSpeciesB = basisCell[0].molecules.getMolecule(molIndex[nCellBox][randomMol1]);
        			pairAB.atom1 = molSpeciesB;
        			uCorrect += potential.energy(pairAB)/2;
        		}
        		
        	}
        	
        }
        
        energyMeter.setTarget(molecule1);
        uOld += energyMeter.getDataAsScalar();
        uOld += uCorrect;
        pair.atom0 = molecule0;
        pair.atom1 = molecule1;
        uOld -= potential.energy(pair);
        
        
        if(uOld > 1e10){
            throw new ConfigurationOverlapException(box);
        }
        
        groupTransVect.setRandomCube(random);
        groupTransVect.TE(stepSize);
        for(int i=0; i<molIndex.length; i++){
        	moveMoleculeAction.actionPerformed(basisCell[0].molecules.getMolecule(molIndex[i][randomMol0]));
        }
        groupTransVect.TE(-1.0);
        for(int i=0; i<molIndex.length; i++){
        	moveMoleculeAction.actionPerformed(basisCell[0].molecules.getMolecule(molIndex[i][randomMol1]));
        }
  
        uCorrect = 0;
        for(int i=0; i<molIndex[0].length; i++){
        	molSpeciesA = basisCell[0].molecules.getMolecule(molIndex[13][i]);
        	pairAB.atom0 = molSpeciesA;
        	
        	if(molSpeciesA == molecule0 || molSpeciesA == molecule1) continue;
        	
        	for(int nCellBox=0; nCellBox<molIndex.length; nCellBox++){
        		if(nCellBox != 13){
        			molSpeciesB = basisCell[0].molecules.getMolecule(molIndex[nCellBox][randomMol0]);
        			pairAB.atom1 = molSpeciesB;
        			uCorrect += potential.energy(pairAB)/2;
        			
        			molSpeciesB = basisCell[0].molecules.getMolecule(molIndex[nCellBox][randomMol1]);
        			pairAB.atom1 = molSpeciesB;
        			uCorrect += potential.energy(pairAB)/2;
        		}
        		
        	}
        	
        }
        
        
        energyMeter.setTarget(molecule0);
        uNew = energyMeter.getDataAsScalar();
        uNew += uCorrect;
        energyMeter.setTarget(molecule1);
        uNew += energyMeter.getDataAsScalar();
        uNew -= potential.energy(pair);
        /*
         * Because we have uNew is infinity, and we don't want to have to 
         * worry about the system subtracting infinity from infinity, and
         * setting uNew equal to zero, and accepting the move.
         */
       // if(Double.isInfinite(uNew)) {return true;}  
       
        
        return true;
    }

    public double getA() {
        return 1.0;
    }

    public double getB() {
        return -(uNew - uOld);
    }

    public void rejectNotify() {
    	
        for(int i=0; i<molIndex.length; i++){
        	moveMoleculeAction.actionPerformed(basisCell[0].molecules.getMolecule(molIndex[i][randomMol0]));
        }
        groupTransVect.TE(-1.0);

        for(int i=0; i<molIndex.length; i++){
        	moveMoleculeAction.actionPerformed(basisCell[0].molecules.getMolecule(molIndex[i][randomMol1]));
        }
    }
    
    private static final long serialVersionUID = 1L;
    protected final MoleculeChildAtomAction moveMoleculeAction;
    protected final IVectorRandom groupTransVect;
    protected IMolecule molecule0, molecule1;
    protected final MeterPotentialEnergy energyMeter;
    protected MoleculeSource moleculeSource;
    protected double uOld, uNew;
    protected final IRandom random;
    protected final AtomIteratorArrayListSimple affectedMoleculeIterator;
    protected final AtomArrayList affectedMoleculeList;
    protected final AtomActionTranslateBy singleAction;
    protected final MoleculePair pair, pairAB;
    protected IPotentialMolecular potential;
    protected int[][] molIndex;
    protected int randomMol0, randomMol1;
    protected IMolecule molSpeciesA, molSpeciesB;
    protected int nA;
    protected BasisCell[] basisCell;
    
}