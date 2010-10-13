package etomica.models.nitrogen;

import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IPotentialMaster;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.space.ISpace;
import etomica.units.Degree;
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
    protected boolean doScaling = true;
	private IVectorMutable[][] initMolecOrientation;
	private IVectorMutable molecOrientation;
    
    public MeterTargetRPMolecule(IPotentialMaster potentialMasterSampled, IPotentialMaster[] potentialMasterMeasured, ISpecies species, ISpace space, ISimulation sim, CoordinateDefinitionNitrogen coordinateDef) {
        this.potentialMasterSampled = potentialMasterSampled;
        this.potentialMasterMeasured = potentialMasterMeasured;
        this.coordinateDefinition = coordinateDef;
        
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
        
        int numMolec = sim.getBox(0).getNMolecules(species);
        molecOrientation = space.makeVector();
    	initMolecOrientation = new IVectorMutable[numMolec][3];
    	/*
		 * initializing the initial orientation of the molecule
		 */
		for (int i=0; i<numMolec; i++){
			initMolecOrientation[i] = space.makeVectorArray(3);
			initMolecOrientation[i] = coordinateDefinition.getMoleculeOrientation(sim.getBox(0).getMoleculeList().getMolecule(i));
		}
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
        
        	if(doScaling){
        		double fac = Math.sqrt((1-Math.cos(Degree.UNIT.toSim(otherAngles[i])))
  						/(1-Math.cos(Degree.UNIT.toSim(angle))));;
        		/*
                 * Re-scaling the coordinate deviation
                 */
          		for (int iMolecule=0; iMolecule<molecules.getMoleculeCount(); iMolecule++){
          			//NOT scaling the translational DOF
          			for(int k=0; k<3; k++){
          				newU[iMolecule*5+k] = u[iMolecule*5+k];
          			}
          			
          			// Scaling the rotational angle for the beta-phase
          			double check = u[(iMolecule*5)+3]*u[(iMolecule*5)+3] + u[(iMolecule*5)+4]*u[(iMolecule*5)+4];
          			
          			//if theta is less or equal to 90deg 
          			// 2.0 is from eq: u3^2 + u4^2 = 2*(1 -cos(theta)
          			
          			//System.out.print("check: "+ check);
          			if(check <= 2.0){
          	        	//System.out.println(" less 2.0");
          				newU[(iMolecule*5)+3] = fac*u[(iMolecule*5)+3];
          				newU[(iMolecule*5)+4] = fac*u[(iMolecule*5)+4];
          				
          			} else if(check>2.0 && check <=4.0) {
          			
          				IMolecule molecule = molecules.getMolecule(iMolecule);
          				IVectorMutable leafPos0 = molecule.getChildList().getAtom(0).getPosition();
          				IVectorMutable leaftPos1 = molecule.getChildList().getAtom(1).getPosition();
          				//Flipping the molecule
          				molecOrientation.Ev1Mv2(leaftPos1, leafPos0);
          				molecOrientation.TE(-1);
          				molecOrientation.normalize();
          				
          				double [] uCalc = calcUMethod(molecOrientation, iMolecule);
          				
          				//System.out.println(" ***********greater 2.0");
          				newU[(iMolecule*5)+3] = fac*uCalc[0];
          				newU[(iMolecule*5)+4] = fac*uCalc[1];
          				
          			} 
          			
          			if((newU[(iMolecule*5)+3]*newU[(iMolecule*5)+3] + newU[(iMolecule*5)+4]*newU[(iMolecule*5)+4]) >4.0){
          				System.out.println("<MeterTargetRPMolecule> newU3^2+newU4^2 is GREATER THAN 4.0");
          				System.out.println("otherAngle: " + otherAngles[i]);
          				System.out.println("[newU3^2+newU4^2]: "+(newU[(iMolecule*5)+3]*newU[(iMolecule*5)+3] + newU[(iMolecule*5)+4]*newU[(iMolecule*5)+4]));
          				System.out.println("u[3 & 4]: "+u[(iMolecule*5)+3]+" "+u[(iMolecule*5)+4]);
          				System.out.println("newU[3 & 4]: "+newU[(iMolecule*5)+3]+" "+newU[(iMolecule*5)+4]);
          				throw new RuntimeException("SO BUSTED!!");
              		}
          		}
        	
        	
        	} else {
         		for (int iCoord=0; iCoord<coordinateDefinition.getCoordinateDim(); iCoord++){
        			newU[iCoord] = u[iCoord];
        	
          		}
        	}
            double otherEnergy = 0;
     
            coordinateDefinition.setToU(pretendMolecules, newU);
            meterPotentialMeasured[i].setBox(pretendBox);
            otherEnergy = meterPotentialMeasured[i].getDataAsScalar();
            //System.out.println((energy-latticeEnergy)+" "+(otherEnergy-latticeEnergy) + " " + (otherEnergy-energy));
            
            double ai = (otherEnergy-latticeEnergy)/temperature;
            //System.out.println(fac+" "+a0+ " " + ai + " "+ (ai-a0));
            
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

    private double[] calcUMethod(IVector vector, int index){
    	
    	double[] u = new double[2];
    	double u3 = vector.dot(initMolecOrientation[index][1]);
		double u4 = vector.dot(initMolecOrientation[index][2]);
		double ratio = Math.abs(u3/u4);
		
		double a = vector.dot(initMolecOrientation[index][0]);
		double theta = Math.acos(a);
		
		if(Math.abs(u4) > -1e-10 && Math.abs(u4) < 1e-10){
			u[0] = Math.sqrt(2*(1-Math.cos(theta)));
			if(u3 <0.0){
				u[0] = -u[0];
			}
			u[1] = u4;
		} else {
			if(u4 < 0.0){
				u[1] = -Math.sqrt(2*(1-Math.cos(theta))/(ratio*ratio+1));
			} else {
				u[1] = Math.sqrt(2*(1-Math.cos(theta))/(ratio*ratio+1));
			}

			if (u3 < 0.0){
				u[0] = -ratio*Math.sqrt(2*(1-Math.cos(theta))/(ratio*ratio+1));
			} else {
				u[0] = ratio*Math.sqrt(2*(1-Math.cos(theta))/(ratio*ratio+1));
			}
		}
		return u; 
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
    
    public boolean isDoScaling() {
		return doScaling;
	}

	public void setDoScaling(boolean doScaling) {
		this.doScaling = doScaling;
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
