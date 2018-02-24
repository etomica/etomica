/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicBcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveTetragonal;
import etomica.normalmode.BasisBigCell;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMolecular;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;



/**
 * 
 * Determine c/a ratio for tetragonal unit cell that minimizes 
 *  the lattice energy of gamma-N2 crystal structure
 * 
 * gamma-N2 crystal Structure
 * 
 * @author Tai Boon Tan
 *
 */
public class MinimizeGammaNitrogenLatticeParameter extends Simulation{

	public MinimizeGammaNitrogenLatticeParameter(Space space, int numMolecule, double density, double ratio) {
        super(space);
        this.space = space;
        this.density = density;
        this.numMolecule = numMolecule;

        nCell = (int) Math.round(Math.pow((numMolecule / 2), 1.0 / 3.0));

        this.a = Math.pow(numMolecule / (ratio * density), 1.0 / 3.0) / nCell;
        this.c = ratio * a;

        potentialMaster = new PotentialMaster();

        Basis basisBCC = new BasisCubicBcc();
        Basis basis = new BasisBigCell(space, basisBCC, new int[]{nCell, nCell, nCell});

        species = new SpeciesN2ShellModel(space);
        addSpecies(species);

        Boundary boundary = new BoundaryRectangularPeriodic(space, new double[]{nCell * a, nCell * a, nCell * c});
        box = this.makeBox(boundary);
        box.setNMolecules(species, numMolecule);
        int[] nCells = new int[]{1, 1, 1};
        Primitive primitive = new PrimitiveTetragonal(space, nCell * a, nCell * c);

        coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
        coordinateDef.setIsGamma();
        coordinateDef.setOrientationVectorGamma(space);
        coordinateDef.initializeCoordinates(nCells);
        double rC = box.getBoundary().getBoxSize().getX(0) * 0.485;
        //System.out.println("Truncation Radius: " + rC);
        potential = new P2NitrogenShellModel(space, rC);
        potential.setBox(box);

        potentialMaster.addPotential(potential, new ISpecies[]{species, species});


    }
	
	public void setRatio(double ratio, double density){
		a = Math.pow(numMolecule/(ratio*density), 1.0/3.0)/nCell;
		c = ratio*a;
		
		int [] nCells = new int[]{1,1,1};
		
		Boundary boundary = new BoundaryRectangularPeriodic(space, new double[]{nCell*a, nCell*a, nCell*c});
		Primitive primitive = new PrimitiveTetragonal(space, nCell*a, nCell*c);
		box.setBoundary(boundary);
		
		//System.out.println("a: " + a + " ; c: " + c);
		Basis basisBCC = new BasisCubicBcc();
		Basis basis = new BasisBigCell(space, basisBCC, new int[]{nCell, nCell, nCell});
		
		coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
		coordinateDef.setIsGamma();
		coordinateDef.setOrientationVectorGamma(space);
		coordinateDef.initializeCoordinates(nCells);
		
		
	}
	
public double findOptRatio(double minRatio, double maxRatio){
    	
    	int bootstrap = 0;
    	double[] u = new double[3];
    	double[] allRatio = new double[3];
       	
    	double ratio = minRatio;
    	
    	MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster, box);
    	
        while (true) {
        	
            setRatio(ratio, density);

            latticeEnergy = meterPE.getDataAsScalar();
            //System.out.println(ratio + " ;a: " + a+ " ;c: " + c +  " ;lattice energy: " + latticeEnergy/numMolecule);
      
            if (bootstrap < 3) {
                allRatio[bootstrap] = ratio;
                u[bootstrap] = latticeEnergy;
                bootstrap++;
               
                ratio += 0.5*(maxRatio-minRatio);
            }
            else {
                if (ratio > allRatio[2]) {
                    allRatio[0] = allRatio[1];
                    allRatio[1] = allRatio[2];
                    allRatio[2] = ratio;
                    u[0] = u[1];
                    u[1] = u[2];
                    u[2] = latticeEnergy;
                }
                else if (ratio < allRatio[0]) {
                    allRatio[2] = allRatio[1];
                    allRatio[1] = allRatio[0];
                    allRatio[0] = ratio;
                    u[2] = u[1];
                    u[1] = u[0];
                    u[0] = latticeEnergy;
                }
                else if (u[2] > u[0]) {
                    maxRatio = allRatio[2];
                    if (ratio > allRatio[1]) {
                        u[2] = latticeEnergy;
                        allRatio[2] = ratio;
                    }
                    else {
                        u[2] = u[1];
                        allRatio[2] = allRatio[1];
                        u[1] = latticeEnergy;
                        allRatio[1] = ratio;
                    }
                }
                else {
                    minRatio = allRatio[0];
                    if (ratio < allRatio[1]) {
                        u[0] = latticeEnergy;
                        allRatio[0] = ratio;
                    }
                    else {
                        u[0] = u[1];
                        allRatio[0] = allRatio[1];
                        u[1] = latticeEnergy;
                        allRatio[1] = ratio;
                    }
                }

                if (u[1] > u[0] && u[1] > u[2]) {
                    // we found a maximum, due to numerical precision failure
                    // just bail and pretend that the middle point is the global minimum
                	
                	a = Math.pow(numMolecule/(allRatio[1]*density), 1.0/3.0)/nCell;
            		c = allRatio[1]*a;
                	//System.out.println("***"+allRatio[1] + " ;a: " + a+ " ;c: " + c);
                     
                    return allRatio[1];
                }
            }

            if (bootstrap == 3) {
                // now estimate minimum in U from the three points.
                double dc01 = allRatio[1]-allRatio[0];
                double dc12 = allRatio[2]-allRatio[1];
                double du01 = u[1]-u[0];
                double du12 = u[2]-u[1];
                double dudc01 = du01/dc01;
                double dudc12 = du12/dc12;
                double m = (dudc12-dudc01)/(0.5*(dc01+dc12));
                ratio = 0.9*(0.5*(allRatio[1]+allRatio[2]) - dudc12/m) + 0.1*(0.5*(allRatio[0]+allRatio[2]));
                if (ratio == allRatio[1] || ratio == allRatio[2]) {
                    ratio = 0.5*(allRatio[1] + allRatio[2]);
                }
                if (ratio == allRatio[0] || ratio == allRatio[1]) {
                    ratio = 0.5*(allRatio[1] + allRatio[0]);
                }
                if (ratio < minRatio) {
                    ratio = 0.5*(minRatio + allRatio[0]);
                }
                if (ratio > maxRatio) {
                    ratio = 0.5*(maxRatio + allRatio[2]);
                }
                        
                if (ratio == allRatio[0] || ratio == allRatio[1] || ratio == allRatio[2]) {
                    // we converged ratio to numerical precision.
                    //System.out.println("ratio "+ ratio);
                	a = Math.pow(numMolecule/(ratio*density), 1.0/3.0)/nCell;
            		c = ratio*a;
                	//System.out.println("***"+ratio + " ;a: " + a+ " ;c: " + c);
                	
                    return ratio;
                }
            }
        }
    }
	
	public static void main (String[] args){
		
		double a = 3.957;
		double c = 5.109;
		double ratio = c/a;
		int n = 14;
		int numMolecule = n*n*n*2;

		int nUnitCell = (int)Math.round(Math.pow((numMolecule/2), 1.0/3.0));
		double density = numMolecule/(nUnitCell*nUnitCell*nUnitCell*a*a*c);
		
		System.out.println("Determine the lattice parameter for gamma-N2 " +
				"by minimizing the lattice energy for " + numMolecule + " molecules.");
		System.out.println("density: " + density);
		System.out.println("intial a: " + a+ " ;c: " + c);
		System.out.println();
		
		MinimizeGammaNitrogenLatticeParameter minimizer = 
			new MinimizeGammaNitrogenLatticeParameter(Space3D.getInstance(3), numMolecule, density, ratio);
		
		double optRatio = minimizer.findOptRatio(1.20, 1.45);
		System.out.println("ratio: " + optRatio);
		System.out.println("a: " + minimizer.getA());
		System.out.println("c: " + minimizer.getC());
		System.out.println("lattice: " + minimizer.getLatticeEnergy());
		
	}

	public double getA() {
		return a;
	}

	public double getC() {
		return c;
	}

	public double getLatticeEnergy() {
		/*
		 * return energy per particle
		 */
		return latticeEnergy/numMolecule;
	}

	protected double a, c, latticeEnergy;
	protected Box box;
	protected Space space;
	protected int numMolecule;
	protected double density;
	protected int nCell;
	protected PotentialMaster potentialMaster;
	protected PotentialMolecular potential;
	protected CoordinateDefinitionNitrogen coordinateDef;
	protected SpeciesN2ShellModel species;
	private static final long serialVersionUID = 1L;
}
