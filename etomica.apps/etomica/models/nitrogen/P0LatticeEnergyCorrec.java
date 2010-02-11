package etomica.models.nitrogen;

import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.potential.PotentialMolecularSoft;
import etomica.space.ISpace;
import etomica.space.Tensor;
import etomica.units.Kelvin;

/**
 *  The lattice energy correction for Nitrogen model
 *  zero-body interaction
 *  
 *  All the parameters are determined from the fitted equation. 
 *   (uLattice - uLattice_infinite) as function of rho
 * 
 * @author Tai Boon Tan
 *
 */
public class P0LatticeEnergyCorrec implements PotentialMolecularSoft{

	public P0LatticeEnergyCorrec(ISpace space){
		coeff = new double[3];
	}
	
	public double energy(IMoleculeList atoms) {
		double rho = numMolec/box.getBoundary().volume();
		return uCorrection(rho);
	}
	
	public void setBox(IBox p){
		this.box = p;
	   if(species==null){
       	throw new RuntimeException("<MCMoveVolumeN2.java> Must set Species First");
       }
                      
       numMolec = p.getNMolecules(species);
       
       if (numMolec == 32){
       	caseNumMolec = 1;
       	coeff[0] =  0.28814;
       	coeff[1] = -8.04733;
       	coeff[2] = - 489961;
       	
       } else if (numMolec == 108){
       	caseNumMolec = 2;
       	coeff[0] = -1.14599;
       	coeff[1] =  226.707;
       	coeff[2] = - 153552;
       	
       } else if (numMolec == 256){
       	caseNumMolec = 3;
       	coeff[0] = -0.210275;
       	coeff[1] =   34.4034;
       	coeff[2] = - 74007.7;
           
       }else if (numMolec == 500){
       	caseNumMolec = 4;
       	coeff[0] = -0.0760888;
       	coeff[1] =    15.1728;
       	coeff[2] = -  30991.5;
           
       } 
	}
	
    private double uCorrection(double rho){
    	double rho2 = rho * rho;
    	
    	if (caseNumMolec==0){
    		return  0.0;
    	
    	} else {
    		/*
    		 * return the correction energy for the total system
    		 * NOT the correction energy per molecule
    		 * The coeff was fitted with energy in K
    		 * so we have to convert the unit to simulation unit 
    		 */
    		return  Kelvin.UNIT.toSim(numMolec*(coeff[0] + coeff[1]*rho + coeff[2]*rho2));
    	
    	}
    }
	
	public ISpecies getSpecies() {
		return species;
	}

	public void setSpecies(ISpecies species) {
		this.species = species;
	}
	
	public double getRange() {
		return 0;
	}

	public IVector[] gradient(IMoleculeList atoms) {
		return null;
	}

	public IVector[] gradient(IMoleculeList atoms, Tensor pressureTensor) {
		return null;
	}

	public double virial(IMoleculeList atoms) {
		return 0;
	}

	public int nBody() {
		return 0;
	}

	private ISpecies species;
	private IBox box;
	private double[] coeff;
    private int caseNumMolec = 0;
    protected int numMolec;
	private static final long serialVersionUID = 1L;

}
