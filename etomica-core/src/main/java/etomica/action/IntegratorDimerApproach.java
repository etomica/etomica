/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorBox;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.PotentialMaster;
import etomica.space.Vector;

/*
 * Integrator for simulation DimerApproach
 * 
 * Uses the site-site models of Rowley et al (2006) to reproduce potential-energy plots for dimers of 
 * ethanol or methanol molecules along particular approach routes
 * 
 * K.R. Schadel 2008 
 */

public class IntegratorDimerApproach extends IntegratorBox {
	
	public IntegratorDimerApproach(PotentialMaster potentialMaster, Box box) {
		
		super(potentialMaster, 0, box);
		
		atomActionTranslateBy = new AtomActionTranslateBy(this.space);
		atomActionRotateBy = new AtomActionRotateBy(this.space);
		
		//The vectors are needed each step; might as well make them just once
		newOriginB = this.space.makeVector();
        translationVector = this.space.makeVector();
		
	}

	public void initializeCoordinates() {
				
		/* *****************************************************************
		 * *****************************************************************
		 * Initialize coordinates of both molecules in standard orientation 
		 * *****************************************************************
		 * *****************************************************************
		 */
		
        monomerA.getType().initializeConformation(monomerA);
        
        monomerB.getType().initializeConformation(monomerB);
        
        /* *****************************************************************
		 * *****************************************************************
		 * Adjust orientation of monomer B from standard orientation
		 * *****************************************************************
		 * *****************************************************************
		 */
        
        double yaw;
        if (routeParams[0].length == 1) {
        	yaw = routeParams[0][0]*Math.PI/180;
        	
        } else if (routeParams[0].length == 2) {
        	double a = routeParams[0][0];
        	double b = routeParams[0][1];
        	yaw = ( a + b*r ) * Math.PI/180;
        } else {
        	throw new IllegalArgumentException("Incorrect number of elements in yaw array");
        }
        
		double pitch =   routeParams[1][0]*Math.PI/180;
		double roll  =   routeParams[2][0]*Math.PI/180;
	
        atomActionRotateBy.setAngles(roll,pitch,yaw);
        atomGroupActionRotateBy = new MoleculeChildAtomAction(atomActionRotateBy);
		atomGroupActionRotateBy.actionPerformed(monomerB);
		
		printKeyAtomPositions();
		
        /*System.out.println();
        System.out.println("After adjusting orientation: ");
        printKeyAtomPositions();*/
        
        /* *****************************************************************
		 * *****************************************************************
		 * Translate origin of monomer B from global origin
		 * *****************************************************************
		 * *****************************************************************
		 */
		double r0 = r;
        delta = r0/75;
		translationVector();
		translateMonomerB(); 
		
		printKeyAtomPositions();
		
		if (route == 15) {
			checkOHBondLength();
		}
		
		/* *****************************************************************
		 * *****************************************************************
		 * set translation Vector that will be used during each call of doStepInternal
		 * *****************************************************************
		 * *****************************************************************
		 */
		 
		translationVector();
	
	}

    protected void doStepInternal() {

        // inverse sine function used to compute phi becomes undefined at small r
        // Using this if statement causes last few steps to affect no change,
		// but prevents one from having to calculate precise number of steps
		// to avoid NaNs
		
      //  if (r > 0.0) {
        	
        	translateMonomerB();
        	
        	//printKeyAtomPositions();

        	
        	if (route == 15) {
        		checkOHBondLength();
        	}
        	
       // }
		
	}
	
	public void translationVector() {	
		
		/* *****************************************************************
		 * *****************************************************************
		 * Adjust position of second monomer
		 * 
		 * During first call, this moves the monomer from the origin
		 * During subsequent calls, this moves the monomer closer to the origin
		 * 
		 * All atoms are translated by the same vector as the alpha carbon
		 * *****************************************************************
		 * *****************************************************************
		 */
		
		// System.out.println(r);
		
		double theta;
		double phi;
		
		if (routeParams[3].length == 1) {
			
			theta = routeParams[3][0]*Math.PI/180;
			
		} else if (routeParams[3].length == 3) {
			
			double a = routeParams[3][0];
			double b = routeParams[3][1];
			double c = routeParams[3][2];
			
			theta = Math.asin( (a/(r*r)) + (b/r) + c ); 
			
		} else {
			
			throw new IllegalArgumentException("Incorrect number of elements in theta array");
			
		}
		
		if (routeParams[4].length == 1) {
			
			phi = routeParams[4][0]*Math.PI/180;
			
		} else if (routeParams[4].length == 3){
			
			double a = routeParams[4][0];
			double b = routeParams[4][1];
			double c = routeParams[4][2];
			
			phi = Math.asin( (a/(r*r)) + (b/r) + c ); // already in radians
			
		} else {
			
			throw new IllegalArgumentException("Incorrect number of elements in phi array");
			
		}
		
		// System.out.println("phi' = " + phi);
		
		// Rowley et al definition of phi (above) a little different...
		phi = 0.5*Math.PI - phi;
		//phi = -phi;
		
		// System.out.println("phi = " + phi);
		
		// Cartesian coordinates of new alpha-carbon position
		double x = r*Math.cos(theta)*Math.sin(phi);
		double y = r*Math.sin(theta)*Math.sin(phi);
		double z = r*Math.cos(phi);
		
		newOriginB.E( new double [] {x, y, z});
		
		translationVector.Ev1Mv2(newOriginB, atom_aC_B.getPosition());
	
	}
	
	public void translateMonomerB() {
		
		atomActionTranslateBy.setTranslationVector(translationVector);
		atomGroupActionTranslateBy = new MoleculeChildAtomAction(atomActionTranslateBy);
        atomGroupActionTranslateBy.actionPerformed(monomerB);
        
        // Update r for next call 
		r = r-delta;
	
	}
	
	public void setRoute(int route) {
		this.route = route;
	}
	
	public void setRouteParams(double[][] routeParams) {
		this.routeParams = routeParams;
	}
	
	public void setMolecules() {
		moleculeList = box.getMoleculeList();
		monomerA = moleculeList.get(0);
		monomerB = moleculeList.get(1);
	}
	
	public void setImportantAtoms() {
		atomSetA = monomerA.getChildList();
	    atomSetB = monomerB.getChildList();
	    
	    atom_O_A  = atomSetA.get(0);
	    atom_aC_A = atomSetA.get(1);
	    atom_aH_A = atomSetA.get(2);
	    atom_H1_A = atomSetA.get(4);
	    
	    atom_O_B  = atomSetB.get(0);
	    atom_aC_B = atomSetB.get(1);
	    atom_aH_B = atomSetB.get(2);
	    
	}
	
	public IAtom getAtom_O_A()  { return atom_O_A ;}
	public IAtom getAtom_aC_A() { return atom_aC_A;}
	public IAtom getAtom_aH_A() { return atom_aH_A;}
	public IAtom getAtom_H1_A() { return atom_H1_A;}
	
	public IAtom getAtom_O_B()  { return atom_O_B ;}
	public IAtom getAtom_aC_B() { return atom_aC_B;}
	public IAtom getAtom_aH_B() { return atom_aH_B;}
	
	public void printKeyAtomPositions() {
		System.out.println();
	    System.out.println("Position of second monomer's oxygen: " + atom_O_B.getPosition()); 
	    System.out.println("Position of second monomer's alpha hydrogen: " + atom_aH_B.getPosition()); 
		System.out.println("Position of second monomer's alpha carbon: " + atom_aC_B.getPosition());  
	    System.out.println();
	}
	
	public void checkOHBondLength() {
		Vector work = space.makeVector();
        
        work.Ev1Mv2(atom_aC_A.getPosition(), atom_aC_B.getPosition());
        double raCaC = Math.sqrt(work.squared());
        
        System.out.println("aC-aC distance" + raCaC);
        
        work.Ev1Mv2(atom_O_B.getPosition(), atom_O_A.getPosition());
        double rOO = Math.sqrt(work.squared());
        
        System.out.println("O-O distance" + rOO);
        
        work.Ev1Mv2(atom_O_B.getPosition(), atom_aH_A.getPosition());
        double rOH1 = Math.sqrt(work.squared());
        
        double bondOHA = rOO - rOH1;
        
        work.Ev1Mv2(atom_O_A.getPosition(), atom_aH_B.getPosition());
        double rOH2 = Math.sqrt(work.squared());
        double bondOHB = rOO - rOH2;
        
        /*System.out.println();
        System.out.println("bondOH (version A)" + bondOHA);
        System.out.println("bondOH (version B)" + bondOHB);
        System.out.println();*/
	}
	
	private static final long serialVersionUID = 1L;
    
    protected AtomActionTranslateBy atomActionTranslateBy;
    protected final AtomActionRotateBy atomActionRotateBy;
   
    protected MoleculeChildAtomAction atomGroupActionRotateBy;
    protected MoleculeChildAtomAction atomGroupActionTranslateBy;
   
    // vectors used in translateMonomerB()
    protected Vector newOriginB;
    protected Vector translationVector;
    
    protected IMoleculeList moleculeList; // List of monomers in box.
    protected IAtomList atomSetA; // List of sites in monomer A.
    protected IAtomList atomSetB; // List of sites in monomer B
    
	protected IMolecule monomerA;
	protected IMolecule monomerB;
	
	// ID's of sites in monomer A
	protected static IAtom atom_O_A;
    protected static IAtom atom_aC_A;
    protected static IAtom atom_aH_A;
    protected static IAtom atom_H1_A; 
    
    // ID's of sites in monomer B
    protected static IAtom atom_O_B;
    protected static IAtom atom_aC_B;
    protected static IAtom atom_aH_B;
    
    // The distance between the alpha carbons of the two monomers
    protected double r = 16.0;  
    
    // The size of the change in r for each translation of monomer B
	protected double delta;
	
	// The number of the approach route
	protected int route;

	// Route-specific parameters
	protected double[][] routeParams;

}
