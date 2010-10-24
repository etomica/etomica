package etomica.models.nitrogen;

import java.io.FileWriter;
import java.io.IOException;

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
import etomica.nbr.list.molecule.PotentialMasterListMolecular;
import etomica.space.ISpace;
import etomica.units.Kelvin;
import etomica.units.Null;

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
public class MeterTargetTPMolecule implements IEtomicaDataSource {

    protected final MeterPotentialEnergy meterPotential;
    protected final IPotentialMaster potentialMaster;
    protected double latticeEnergy;
    protected double temperature;
    protected double[] otherTemperatures;
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
    protected FileWriter fw;
    protected boolean isBetaPhase = false;
	private IVectorMutable[][] initMolecOrientation;
	private IVectorMutable molecOrientation;
   
    public MeterTargetTPMolecule(IPotentialMaster potentialMaster, ISpecies species, ISpace space, ISimulation sim) {
        this.potentialMaster = potentialMaster;
        meterPotential = new MeterPotentialEnergy(potentialMaster);
        this.species = species;
        pretendBox = new Box(space);
        sim.addBox(pretendBox);
     
        tag = new DataTag();
        
        molecOrientation = space.makeVector();
        
        int numMolec = sim.getBox(0).getNMolecules(species);
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
           		// The alpha-phase Nitrogen
           		
        		/*
                 * Re-scaling the coordinate deviation
                 */
          		for (int iMolecule=0; iMolecule<molecules.getMoleculeCount(); iMolecule++){
          			//NOT scaling the translational DOF
          			for(int k=0; k<3; k++){
          				newU[iMolecule*5+k] = fac*u[iMolecule*5+k];
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
          				System.out.println("<MeterTargetTPMolecule> newU3^2+newU4^2 is GREATER THAN 4.0");
          				System.out.println("[newU3^2+newU4^2]: "+(newU[(iMolecule*5)+3]*newU[(iMolecule*5)+3] + newU[(iMolecule*5)+4]*newU[(iMolecule*5)+4]));
          				System.out.println("u[3 & 4]: "+u[(iMolecule*5)+3]+" "+u[(iMolecule*5)+4]);
          				System.out.println("newU[3 & 4]: "+newU[(iMolecule*5)+3]+" "+newU[(iMolecule*5)+4]);
          				throw new RuntimeException("SO BUSTED!!");
              		}
          		}
        		
           	}
            
            coordinateDefinition.setToU(pretendMolecules, newU);
            otherEnergy = meterPotential.getDataAsScalar();
                        
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
        IBox realBox = coordinateDefinition.getBox();
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
}
