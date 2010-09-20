package etomica.models.nitrogen;

import java.io.FileWriter;
import java.io.IOException;

import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.api.IPotentialMaster;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.space.ISpace;
import etomica.units.Kelvin;
import etomica.units.Null;

/**
 * Rotational Perturbation for beta-phase Nitrogen
 * - scaling of the rotational angle to better map the phase space
 * of the perturbing systems.
 * 
 * 
 * ep/(e0 + alpha*ep)  if  Tp > T0
 * e0/(ep + alpha*e0)  if  T0 > Tp
 * 
 * @author taitan
 *
 */
public class MeterTargetRPMolecule implements IEtomicaDataSource {

    protected MeterPotentialEnergy meterPotentialSampled; 
    protected MeterPotentialEnergy[] meterPotentialMeasured;
    protected IPotentialMaster potentialMasterSampled;
    protected IPotentialMaster[] potentialMasterMeasured;
    protected double latticeEnergy;
    protected double temperature;
    protected double angle;
    protected double[] otherAngles;
    protected DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected final DataTag tag;
    protected final IBox pretendBox;
    protected CoordinateDefinitionNitrogen coordinateDefinition;
    protected final ISpecies species;
    protected double[][] alpha;
    protected double[] alphaCenter;
    protected double alphaSpan;
    protected int numAlpha = 1;
    
    public MeterTargetRPMolecule(IPotentialMaster potentialMasterSampled, IPotentialMaster[] potentialMasterMeasured, ISpecies species, ISpace space, ISimulation sim, CoordinateDefinitionNitrogen coordinateDef) {
        this.potentialMasterSampled = potentialMasterSampled;
        this.potentialMasterMeasured = potentialMasterMeasured;
        this.coordinateDefinition = coordinateDef;
        
        System.out.println("potentialMaster: " + potentialMasterMeasured[0]);
        System.out.println("num: " + potentialMasterMeasured.length);
        meterPotentialSampled = new MeterPotentialEnergy(potentialMasterSampled);
        meterPotentialMeasured = new MeterPotentialEnergy[potentialMasterMeasured.length];
        
        for (int i=0; i<meterPotentialMeasured.length; i++){
        	meterPotentialMeasured[i] = new MeterPotentialEnergy(potentialMasterMeasured[i]);
        }
        this.species = species;
        pretendBox = new Box(space);
        pretendBox.setBoundary(coordinateDef.getBox().getBoundary());
        
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
    	IBox realBox = coordinateDefinition.getBox();
        meterPotentialSampled.setBox(realBox);
        
        double energy = meterPotentialSampled.getDataAsScalar();
    
        
        pretendBox.setBoundary(realBox.getBoundary());
        pretendBox.setNMolecules(species, realBox.getNMolecules(species));
        
        IMoleculeList molecules = realBox.getMoleculeList();
        IMoleculeList pretendMolecules = pretendBox.getMoleculeList();
        
        double a0 = (energy-latticeEnergy)/temperature;
        double[] x = data.getData();
        
        double[] u = coordinateDefinition.calcU(molecules);
        double[] newU = new double[coordinateDefinition.getCoordinateDim()];
        
        for (int i=0; i<otherAngles.length; i++) {
            double fac = (otherAngles[i]/angle);
            double otherEnergy = 0;
            /*
             * Re-scaling the coordinate deviation
             */
      		for (int iCoord=0; iCoord<coordinateDefinition.getCoordinateDim(); iCoord++){
      			// Scaling the rotational angle for the beta-phase
      			if(iCoord>0 && (iCoord%5==3 || iCoord%5==4)){
      				newU[iCoord] = fac*u[iCoord];
      			} else {
      				newU[iCoord] = u[iCoord];
      			}
      		}

            coordinateDefinition.setToU(pretendMolecules, newU);
            meterPotentialMeasured[i].setBox(pretendBox);
            otherEnergy = meterPotentialMeasured[i].getDataAsScalar();
            System.out.println((energy-latticeEnergy)+" "+(otherEnergy-latticeEnergy) + " " + (otherEnergy-energy));
            
            double ai = (otherEnergy-latticeEnergy)/temperature;
           // System.out.println(fac+" "+a0+ " " + ai + " "+ (ai-a0));
            
            for (int j=0; j<numAlpha; j++) {
                if (angle> otherAngles[i]) {
                    x[i*numAlpha+j] = 1.0/(alpha[i][j]+Math.exp(ai-a0));
                }
                else {
                    x[i*numAlpha+j] = 1.0/(1+alpha[i][j]*Math.exp(ai-a0));
                }
            }
        }
        
        return data;
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

    public double[] getOtherAngles() {
        return otherAngles;
    }

    public void setOtherAngles(double[] otherAngles) {
        this.otherAngles = otherAngles;
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

	public double getAngle() {
		return angle;
	}

	public void setAngle(double angle) {
		this.angle = angle;
	}
    

}
