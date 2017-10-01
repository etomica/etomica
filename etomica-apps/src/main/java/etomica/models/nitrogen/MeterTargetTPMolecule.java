/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataSource;
import etomica.data.IEtomicaDataInfo;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.molecule.IMoleculeList;
import etomica.nbr.list.molecule.PotentialMasterListMolecular;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.units.Kelvin;
import etomica.units.dimensions.Null;

import java.io.FileWriter;
import java.io.IOException;

/**
 * For molecular model, e.g., diatomic Nitrogen model, with 5 d.o.f.
 *  or even molecular model with 6 d.o.f.
 *  d.o.f. = degrees of freedom. 
 * 
 * 
 * Meter that measures the overlap averages for perturbing from a solid at one
 * temperature into other the system at other temperatures.  The atoms are
 * scaled back toward their lattice sites by a factor of Tp/T0, where T0 is the
 * simulation temperature and Tp is the temperature being perturbed into.
 * 
 * For the purposes of the overlap average, the lower temperature is considered
 * the reference.
 * 
 * ep/(e0 + alpha*ep)  if  Tp > T0
 * e0/(ep + alpha*e0)  if  T0 > Tp
 * 
 * @author taitan
 *
 */
public class MeterTargetTPMolecule implements IDataSource {

    protected final MeterPotentialEnergy meterPotential;
    protected final PotentialMaster potentialMaster;
    protected double latticeEnergy;
    protected double temperature;
    protected double[] otherTemperatures;
    protected DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected final DataTag tag;
    protected final Box pretendBox;
    protected CoordinateDefinitionNitrogen coordinateDefinition;
    protected final ISpecies species;
    protected double[][] alpha;
    protected double[] alphaCenter;
    protected double alphaSpan;
    protected int numAlpha = 1;
    protected FileWriter fw;
    protected boolean isBetaPhase = false;
   
    public MeterTargetTPMolecule(PotentialMaster potentialMaster, ISpecies species, Space space, Simulation sim) {
        this.potentialMaster = potentialMaster;
        meterPotential = new MeterPotentialEnergy(potentialMaster);
        this.species = species;
        pretendBox = new Box(space);
        sim.addBox(pretendBox);
     
        tag = new DataTag();
    }

	public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    public IData getData() {
    	Box realBox = coordinateDefinition.getBox();
        meterPotential.setBox(realBox);
        double energy = meterPotential.getDataAsScalar();
        meterPotential.setBox(pretendBox);
 
        pretendBox.setBoundary(realBox.getBoundary());      
        IMoleculeList molecules = realBox.getMoleculeList();
        IMoleculeList pretendMolecules = pretendBox.getMoleculeList();
 
        double a0 = (energy-latticeEnergy)/temperature;
        double[] x = data.getData();
        
        double[] u = coordinateDefinition.calcU(molecules);
        double[] newU = new double[coordinateDefinition.getCoordinateDim()];

        for (int i=0; i<otherTemperatures.length; i++) {
            double fac = Math.sqrt(Kelvin.UNIT.toSim(otherTemperatures[i])/temperature);
            double otherEnergy = 0;
            /*
             * Re-scaling the coordinate deviation
             */
          
        	boolean notOverScale = true;
        	
          	if(isBetaPhase){
          		for (int iCoord=0; iCoord<coordinateDefinition.getCoordinateDim(); iCoord++){
          			// NOT Scaling the rotational angle for the beta-phase

          			if(iCoord>0 && (iCoord%5==3 || iCoord%5==4)){
          				newU[iCoord] = u[iCoord];
                    } else {
                    	newU[iCoord] = fac*u[iCoord];
                    }
          		}
           	} else {
           		for (int iCoord=0; iCoord<coordinateDefinition.getCoordinateDim(); iCoord++){
          			newU[iCoord] = fac*u[iCoord];
                }	
           		
           	  	double totalCosTheta = 0.0;
              	
              	for(int iU=0; iU<newU.length; iU+=5){
              		double u3 = newU[iU+3];
              		double u4 = newU[iU+4];
              		
              		double costheta = 1- 0.5*(u3*u3 + u4*u4);
              		
              		totalCosTheta += costheta;
              		
              		// checking for the validity of orientation scaling
              		if(costheta < 0.0 ){
//              			System.out.println("***** "+costheta + " " +u3+" " + u4);
              			otherEnergy = Double.POSITIVE_INFINITY;
              			notOverScale = false;
              			break;
              		}
              	}
              	
              	// checking for the disorderness of alpha phase
              	if(notOverScale){
              		double aveCosTheta = totalCosTheta/molecules.getMoleculeCount();
                  	if(aveCosTheta < 0.8 ){
                  		otherEnergy = Double.POSITIVE_INFINITY;
                  		notOverScale = false;
                  	}  
              	}
              	
           	}
            
          	if(notOverScale){
          		coordinateDefinition.setToU(pretendMolecules, newU);
            	otherEnergy = meterPotential.getDataAsScalar();
          	} 
        
            double ai = (otherEnergy-latticeEnergy)/Kelvin.UNIT.toSim(otherTemperatures[i]);
            //System.out.println("ai-a0: " + ai + " " + a0 + " "+ (ai-a0));
            
            for (int j=0; j<numAlpha; j++) {
                if (temperature>Kelvin.UNIT.toSim(otherTemperatures[i])) {
                    x[i*numAlpha+j] = 1.0/(alpha[i][j]+Math.exp(ai-a0));
                }
                else {
                    x[i*numAlpha+j] = 1.0/(1+alpha[i][j]*Math.exp(ai-a0));
                }
            }
        }
        if (fw != null) {
            try {
                fw.write(x[(numAlpha-1)/2]+"\n");
            }
            catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
        return data;
    }


    /**
     * Writes collected overlap data (for the "middle" alpha) to a file.
     * Only data for the first perturbed temperature is written.
     */
    public void openFW(String filename) {
        try {
            if (fw != null) {
                fw.close();
            }
            fw = new FileWriter(filename);
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Closes file with overlap data.
     */
    public void closeFW() {
        try {
            fw.close();
            fw = null;
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public double getLatticeEnergy() {
        return latticeEnergy;
    }

    public void setLatticeEnergy(double latticeEnergy) {
        this.latticeEnergy = latticeEnergy;
    }

    public double getTemperature() {
        return temperature;
    }

    public void setTemperature(double temperature) {
        this.temperature = temperature;
    }

    public double[] getOtherTemperatures() {
        return otherTemperatures;
    }

    public void setOtherTemperatures(double[] otherTemperatures) {
        this.otherTemperatures = otherTemperatures;
    }
    
    protected void initAlpha() {
        if (alphaCenter == null) {
            return;
        }
        alpha = new double[alphaCenter.length][numAlpha];
        for (int i=0; i<alpha.length; i++) {
            if (numAlpha == 1) {
                alpha[i][0] = alphaCenter[i];
            }
            else {
                for (int j=0; j<numAlpha; j++) {
                    alpha[i][j] = alphaCenter[i]*Math.exp(2.0*alphaSpan*(j-(numAlpha-1)/2)/(numAlpha-1));
                }
            }
        }
        data = new DataDoubleArray(numAlpha*alphaCenter.length);
        dataInfo = new DataInfoDoubleArray("overlap", Null.DIMENSION, new int[]{numAlpha*alphaCenter.length});
    }
    
    public double[] getAlpha(int iTemp) {
        return alpha[iTemp];
    }
    
    public void setAlpha(double[] newAlpha) {
        alphaCenter = newAlpha;
        initAlpha();
    }
    
    public void setAlphaSpan(double newAlphaSpan) {
        alphaSpan = newAlphaSpan;
        initAlpha();
    }
    
    public void setNumAlpha(int newNumAlpha) {
        numAlpha = newNumAlpha;
        initAlpha();
    }
    
    
    public boolean isBetaPhase() {
		return isBetaPhase;
	}

	public void setBetaPhase(boolean isBetaPhase) {
		this.isBetaPhase = isBetaPhase;
	}

    public CoordinateDefinitionNitrogen getCoordinateDefinition() {
        return coordinateDefinition;
    }

    public void setCoordinateDefinition(CoordinateDefinitionNitrogen newCoordinateDefinition) {
        this.coordinateDefinition = newCoordinateDefinition;

        // insert molecules into the box at their lattice sites.
        // we do this because want to find neighbors now (and then never again)
        Box realBox = coordinateDefinition.getBox();
        pretendBox.setBoundary(realBox.getBoundary());
        pretendBox.setNMolecules(species, realBox.getNMolecules(species));
        IMoleculeList pretendMolecules = pretendBox.getMoleculeList();
        
        double[] u = new double[coordinateDefinition.getCoordinateDim()];
        coordinateDefinition.setToU(pretendMolecules, u);

        if (potentialMaster instanceof PotentialMasterListMolecular) {
            // find neighbors now.
            ((PotentialMasterListMolecular)potentialMaster).getNeighborManager(pretendBox).reset();
        }
    }

}
