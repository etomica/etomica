package etomica.normalmode;

import etomica.atom.Atom;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.molecule.MoleculeAgentManager;
import etomica.space.Space;
import etomica.space.Vector;
import org.apache.commons.math3.special.Erf;
import etomica.atom.AtomLeafAgentManager;
import etomica.potential.EwaldSummation;

/**
@author  Weisong Lin
 */
public class EwaldSummationLattice extends EwaldSummation {
	protected final Vector ldr;
	protected final Vector shift;
	protected final MoleculeAgentManager latticeCoordinates;
	protected IAtomPositionDefinition positionDefinition;
	protected final Vector rA;
	
	public EwaldSummationLattice(Box box,
			AtomLeafAgentManager<MyCharge> atomAgentManager, Space _space,
			double kCut, double rCutRealES) {
		this(box,atomAgentManager,_space,kCut,rCutRealES,null);
		
	}
	
	public EwaldSummationLattice(Box box,
			AtomLeafAgentManager<MyCharge> atomAgentManager, Space _space,
			double kCut, double rCutRealES,MoleculeAgentManager latticeCoordinates) {
		super(box, atomAgentManager, _space, kCut, rCutRealES);
		
		ldr = space.makeVector();
		shift = space.makeVector();
		rA = space.makeVector();
		this.latticeCoordinates = latticeCoordinates;
		this.positionDefinition = new AtomPositionCOM(space);
	}
	
	  public double uReal(){
		  
	        int nAtoms = box.getLeafList().getAtomCount();
	        double uReal = 0.0;
	        for (int i=0; i < nAtoms; i++){//H
	            Atom atomA = box.getLeafList().getAtom(i);
	            double chargeA = atomAgentManager.getAgent(atomA).charge;
	            rA.E(latticeCoordinates == null?positionDefinition.position(atomA.getParentGroup()):
	            	((MoleculeSiteSource.LatticeCoordinate)latticeCoordinates.getAgent(atomA.getParentGroup())).position);
            	
	            if (chargeA==0) continue;
	 
	            int aIndex = atomA.getParentGroup().getIndex();
	            Vector positionA = atomA.getPosition();
	            for (int j=i; j < nAtoms; j++){
	                Atom atomB = box.getLeafList().getAtom(j);

	                int bIndex = atomB.getParentGroup().getIndex();
	     
	                if(aIndex == bIndex && nRealShells[0]==0 && nRealShells[1]==0 && nRealShells[2]==0) continue;//Skip atom-pairs in the same molecule in the orig. cell.
	      
	                double chargeB = atomAgentManager.getAgent(atomB).charge;
	              
	                if (chargeB==0) continue;
	                
	                Vector rB = latticeCoordinates == null?positionDefinition.position(atomB.getParentGroup()):
		            	((MoleculeSiteSource.LatticeCoordinate)latticeCoordinates.getAgent(atomB.getParentGroup())).position;
	                ldr.Ev1Mv2(rA, rB);
	                shift.Ea1Tv1(-1, ldr);
	                box.getBoundary().nearestImage(ldr);
	                shift.PE(ldr);
	                Vector positionB = atomB.getPosition();
	                rAB.Ev1Mv2(positionA, positionB);// get vector rAB
	                rAB.PE(shift);
	                
	                // minimum image
	                
	            	boolean isSelf = i == j;
	                for(int nx = -nRealShells[0]; nx <= nRealShells[0]; nx++) {
	                    Lxyz.setX(0, nx*boxSize[0]); 
	                    for(int ny = -nRealShells[1]; ny <= nRealShells[1]; ny++) {
	                        Lxyz.setX(1, ny*boxSize[1]);
	                        for(int nz = -nRealShells[2]; nz <= nRealShells[2]; nz++) {
	                        	boolean centerImage = nx*nx+ny*ny+nz*nz == 0;
	                            if (aIndex==bIndex && centerImage) continue;//Skip atom-pairs in the same molecule in the orig. cell & ignores self+centerImage too
	                            Lxyz.setX(2, nz*boxSize[2]);
	                            drTmp.Ev1Pv2(ldr, Lxyz);
	                            double r2 = drTmp.squared();
	                            if(r2 > rCutSquared) continue;
	                            drTmp.Ev1Pv2(rAB, Lxyz);
	                            r2 = drTmp.squared();
	                            double drTmpM = Math.sqrt(r2);
	                            double tmepReal = chargeA * chargeB * Erf.erfc(alpha * drTmpM) / drTmpM;//Don't worry about 1/2 factor;j>i
	                            uReal+= (isSelf ? 0.5 : 1.0)*tmepReal;
	                            
//	                            if(Math.abs((isSelf ? 0.5 : 1.0)*tmepReal )< 0.0001){
//	        						System.out.println(tmepReal);
//	        						throw new RuntimeException();
//	        					}
	                            
//	                            try{
//	                            	fileWriter.write(  i + " " + j + " " + (isSelf ? 0.5 : 1.0)*tmepReal + "\n");
//	                            }
//	                            catch (IOException e){
//	                        		throw new RuntimeException(e);
//	                        	} 
	                        }
	                    }
	                }
	            }// close for all sites in j-th molecule
	        } // close for the outside loop
//	        System.out.println("uReal = " + uReal);
	        return uReal;
	    }
	  protected void realGradient (IAtomList atoms){
	    	int nAtoms = box.getLeafList().getAtomCount();
	    	//Real gradient  //Cross Interaction
	        for (int i=0; i < nAtoms; i++){
	            IAtom atomA = box.getLeafList().getAtom(i);
	            double chargeA = atomAgentManager.getAgent(atomA).charge;
	            if (chargeA==0) continue;
	            int aIndex = atomA.getParentGroup().getIndex(); // molecule a
	            Vector positionA = atomA.getPosition();
	            rA.E(latticeCoordinates == null?positionDefinition.position(atomA.getParentGroup()):
	            	((MoleculeSiteSource.LatticeCoordinate)latticeCoordinates.getAgent(atomA.getParentGroup())).position);
	            for (int j=i+1; j < nAtoms; j++){
	                IAtom atomB = box.getLeafList().getAtom(j);
	                int bIndex = atomB.getParentGroup().getIndex(); // molecule b

	                if(nRealShells[0] == 0 && nRealShells[1] == 0 && nRealShells[2] == 0 && aIndex == bIndex) continue;//Skip same molecules!

	                double chargeB = atomAgentManager.getAgent(atomB).charge;
	                if (chargeB==0) continue;
	                Vector positionB = atomB.getPosition();
	                Vector rB = latticeCoordinates == null?positionDefinition.position(atomB.getParentGroup()):
		            	((MoleculeSiteSource.LatticeCoordinate)latticeCoordinates.getAgent(atomB.getParentGroup())).position;
	                ldr.Ev1Mv2(rA, rB);
                    shift.Ea1Tv1(-1, ldr);
                    box.getBoundary().nearestImage(ldr);
                    shift.PE(ldr);
	                rAB.Ev1Mv2(positionA, positionB); //rAB == rA - rB
	                rAB.PE(shift);
	                for (int nx = -nRealShells[0]; nx <= nRealShells[0]; nx++) {
	                    Lxyz.setX(0, nx*boxSize[0]); 
	                    for (int ny = -nRealShells[1]; ny <= nRealShells[1]; ny++) {
	                        Lxyz.setX(1, ny*boxSize[1]);
	                        for (int nz = -nRealShells[2]; nz <= nRealShells[2]; nz++) {
	                            if (aIndex==bIndex && nx*nx+ny*ny+nz*nz == 0) continue;
	                            Lxyz.setX(2, nz*boxSize[2]);
	                            
	                            
	                           
	                            drTmp.Ev1Pv2(ldr, Lxyz);
	                            double ldr2 = drTmp.squared();
	                            if (ldr2 > rCutSquared) continue; 
	                            
	                            drTmp.Ev1Pv2(rAB, Lxyz);
	                            
	                            double rAB2 = drTmp.squared();
	                            double rABMagnitude = Math.sqrt(rAB2);
	                            double rAB3 = rABMagnitude*rAB2;
	                            double B = Erf.erfc(alpha*rABMagnitude) + 2.0*alpha*rABMagnitude/sqrtPI * Math.exp(-alpha2*rAB2) ;
	                            double realCoeff = - chargeA*chargeB * B / rAB3; // gradU = -F
	                            gradient[i].PEa1Tv1(realCoeff, drTmp);
	                            gradient[j].PEa1Tv1(-realCoeff, drTmp);
	                        }
	                    }
	                }
	            }
	        }
	    }
	  
	  
}
