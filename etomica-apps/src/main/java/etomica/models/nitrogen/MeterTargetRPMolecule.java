/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.api.ISpecies;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.molecule.IMoleculeList;
import etomica.nbr.list.molecule.PotentialMasterListMolecular;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
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

    protected MeterPotentialEnergy meterPotential; 
    protected PotentialMaster potentialMaster;
    protected double latticeEnergy;
    protected double temperature;
    protected double angle;
    protected double[] otherAngles;
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
    protected boolean doScaling = true;
	private Vector[][] initMolecOrientation;
	protected PRotConstraint pRotConstraint;
    
    public MeterTargetRPMolecule(PotentialMaster potentialMasterSampled, ISpecies species, Space space, Simulation sim,
                                 CoordinateDefinitionNitrogen coordinateDef, PRotConstraint pRotConstraint) {
        this.potentialMaster = potentialMasterSampled;
        this.coordinateDefinition = coordinateDef;
        this.pRotConstraint = pRotConstraint;
        
        meterPotential = new MeterPotentialEnergy(potentialMasterSampled);
        this.species = species;
        pretendBox = new Box(space);
        pretendBox.setBoundary(coordinateDef.getBox().getBoundary());
        
        sim.addBox(pretendBox);
        tag = new DataTag();
        
        int numMolec = sim.getBox(0).getNMolecules(species);
    	initMolecOrientation = new Vector[numMolec][3];
    	/*
		 * initializing the initial orientation of the molecule
		 */
		for (int i=0; i<numMolec; i++){
			initMolecOrientation[i] = space.makeVectorArray(3);
			initMolecOrientation[i] = coordinateDefinition.getMoleculeOrientation(sim.getBox(0).getMoleculeList().getMolecule(i));
		}
		
		Box realBox = coordinateDef.getBox();
		pretendBox.setBoundary(realBox.getBoundary());
        pretendBox.setNMolecules(species, realBox.getNMolecules(species));
        
        IMoleculeList pretendMolecule = pretendBox.getMoleculeList();
        double[] initU = new double[coordinateDef.getCoordinateDim()];
        coordinateDef.setToU(pretendMolecule, initU);
        ((PotentialMasterListMolecular)potentialMasterSampled).getNeighborManager(pretendBox).reset();
		
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
    
        pretendBox.setBoundary(realBox.getBoundary());
        pretendBox.setNMolecules(species, realBox.getNMolecules(species));
        
        IMoleculeList molecules = realBox.getMoleculeList();
        IMoleculeList pretendMolecules = pretendBox.getMoleculeList();
        
        double a0 = (energy-latticeEnergy)/temperature;
        double[] x = data.getData();
        
        double[] u = coordinateDefinition.calcU(molecules);
        double[] newU = new double[u.length];
        
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
          			
          			newU[(iMolecule*5)+3] = fac*u[(iMolecule*5)+3];
      				newU[(iMolecule*5)+4] = fac*u[(iMolecule*5)+4];
          			
          			if((newU[(iMolecule*5)+3]*newU[(iMolecule*5)+3] + newU[(iMolecule*5)+4]*newU[(iMolecule*5)+4]) >2.0){
          				System.out.println("<MeterTargetRPMolecule> newU3^2+newU4^2 is GREATER THAN 2.0");
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
            
            pRotConstraint.setSwitch(false);
            meterPotential.setBox(pretendBox);
            otherEnergy = meterPotential.getDataAsScalar();
            pRotConstraint.setSwitch(true);
            
            double ai = (otherEnergy-latticeEnergy)/temperature;
//            System.exit(1);
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
