/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.action.BoxInflate;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.units.dimensions.Dimension;
import etomica.units.Kelvin;
import etomica.units.dimensions.Pressure;
import etomica.util.random.IRandom;

/**
 * Monte Carlo volume-change move for simulations in the NPT ensemble.
 * for Nitrogen molecule.
 * This class is created to take the energy correction into account. 
 *
 * @author Tai Boon Tan
 */
public class MCMoveVolumeN2 extends MCMoveBoxStep {

    public MCMoveVolumeN2(Simulation sim, PotentialMaster potentialMaster,
                          Space _space) {
        this(potentialMaster, sim.getRandom(), _space, 1.0);
    }
    
    /**
     * @param potentialMaster an appropriate PotentialMaster instance for calculating energies
     * @param space the governing space for the simulation
     */
    public MCMoveVolumeN2(PotentialMaster potentialMaster, IRandom random,
                          Space _space, double pressure) {
        super(potentialMaster);
        this.random = random;
        this.D = _space.D();
        inflate = new BoxInflate(_space);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.10);
        setPressure(pressure);
        energyMeter.setIncludeLrc(true);
        affectedAtomIterator = new AtomIteratorLeafAtoms();
        coeff = new double[3];
        rScale = _space.makeVector();
    }
    
    public void setBox(Box p) {
        super.setBox(p);
        energyMeter.setBox(p);
        inflate.setBox(p);
        affectedAtomIterator.setBox(p);
        if(species==null){
        	throw new RuntimeException("<MCMoveVolumeN2.java> Must set Species First");
        }
                       
        numMolec = p.getNMolecules(species);
        double rho = numMolec/p.getBoundary().volume();
        
        if (numMolec == 32){
        	caseNumMolec = 1;
        	coeff[0] =  0.28814;
        	coeff[1] = -8.04733;
        	coeff[2] = - 489961;
        	setLatticeCorrec(uCorrection(rho));
        	
        } else if (numMolec == 108){
        	caseNumMolec = 2;
        	coeff[0] = -1.14599;
        	coeff[1] =  226.707;
        	coeff[2] = - 153552;
        	setLatticeCorrec(uCorrection(rho));
        } else if (numMolec == 256){
        	caseNumMolec = 3;
        	coeff[0] = -0.210275;
        	coeff[1] =   34.4034;
        	coeff[2] = - 74007.7;
         	setLatticeCorrec(uCorrection(rho));
            
        }else if (numMolec == 500){
        	caseNumMolec = 4;
        	coeff[0] = -0.0760888;
        	coeff[1] =    15.1728;
        	coeff[2] = -  30991.5;
        	setLatticeCorrec(uCorrection(rho));
            
        }
        
    }
    
    public boolean doTrial() {
        double vOld = box.getBoundary().volume();
        
        double rhoOld = numMolec/vOld;
        uOld = energyMeter.getDataAsScalar() + uCorrection(rhoOld);
        hOld = uOld + pressure*vOld;
        
        if(isVolChange){
	        vScale = (2.*random.nextDouble()-1.)*stepSize;
	        vNew = vOld * Math.exp(vScale); //Step in ln(V)
	        double scale = Math.exp(vScale/D);
	        rScale.E(new double[]{scale,scale,scale});
	        
        }
        
        if (isXYZChange){
            double xOld = box.getBoundary().getBoxSize().getX(0);
            double yOld = box.getBoundary().getBoxSize().getX(1);
            double zOld = box.getBoundary().getBoxSize().getX(2);
        	xScale = Math.exp((2.*random.nextDouble()-1.)*stepSize);
        	yScale = Math.exp((2.*random.nextDouble()-1.)*stepSize);
        	zScale = Math.exp((2.*random.nextDouble()-1.)*stepSize);
        	
        	vNew = (xOld*xScale)*(yOld*yScale)*(zOld*zScale);
        	rScale.E(new double[]{xScale,yScale,zScale});
        }
        
        rhoNew = numMolec/vNew;
        inflate.setVectorScale(rScale);
        inflate.actionPerformed();
        
        uNew = energyMeter.getDataAsScalar() + uCorrection(rhoNew);
        hNew = uNew + pressure*vNew;
        return true;
    }//end of doTrial
    
    public double getA() {
        return Math.exp((box.getMoleculeList().getMoleculeCount()+1)*vScale);
    }
    
    public double getB() {
        return -(hNew - hOld);
    }
    
    public void acceptNotify() {  /* do nothing */}
    
    public void rejectNotify() {
        inflate.undo();
    }

    public double energyChange() {return uNew - uOld;}
    
    public AtomIterator affectedAtoms() {
        return affectedAtomIterator;
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
	
    public double getLatticeCorrec() {
    	/*
    	 * return latticeCorrec in sim unit
    	 * per molecule
    	 */
		return latticeCorrec;
	}

	public void setLatticeCorrec(double latticeCorrec) {
		this.latticeCorrec = latticeCorrec;
	}
	
	public void setXYZChange(){
		isVolChange = false;
		isXYZChange = true;
	}
	
	private Vector rScale;
	private double[] coeff;
    private int caseNumMolec = 0;
    private ISpecies species;
    protected int numMolec;
	protected double latticeCorrec;
	
    private static final long serialVersionUID = 2L;
    protected double pressure;
    private MeterPotentialEnergy energyMeter;
    protected final BoxInflate inflate;
    private final int D;
    private IRandom random;
    protected final AtomIteratorLeafAtoms affectedAtomIterator;

    private transient double uOld, hOld, vNew, vScale, hNew, rhoNew, xScale, yScale, zScale;
    private transient double uNew = Double.NaN;
    private boolean isVolChange = true;
    private boolean isXYZChange = false;
	
	public void setPressure(double p) {pressure = p;}
    public final double getPressure() {return pressure;}
    public Dimension getPressureDimension() {return Pressure.DIMENSION;}
}
