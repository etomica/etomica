package etomica.normalmode;

import java.io.FileWriter;
import java.io.IOException;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.api.IVectorMutable;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.nbr.list.PotentialMasterList;
import etomica.space.ISpace;
import etomica.units.Null;

/**
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
 */
public class MeterTargetTP implements IEtomicaDataSource {

    protected final MeterPotentialEnergy meterPotential;
    protected final IPotentialMaster potentialMaster;
    protected double latticeEnergy;
    protected double temperature;
    protected double[] otherTemperatures;
    protected DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected final DataTag tag;
    protected final IBox pretendBox;
    protected CoordinateDefinition coordinateDefinition;
    protected final ISpecies species;
    protected double[][] alpha;
    protected double[] alphaCenter;
    protected double alphaSpan;
    protected int numAlpha = 1;
    public static FileWriter fw;
    
    public MeterTargetTP(IPotentialMaster potentialMaster, ISpecies species, ISpace space, ISimulation sim) {
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
        IBox realBox = coordinateDefinition.getBox();
        meterPotential.setBox(realBox);
        double u = meterPotential.getDataAsScalar();
        meterPotential.setBox(pretendBox);

        pretendBox.setBoundary(realBox.getBoundary());
        IAtomList atoms = realBox.getLeafList();
        IAtomList pretendAtoms = pretendBox.getLeafList();
        double a0 = (u-latticeEnergy)/temperature;

        double[] x = data.getData();
        for (int i=0; i<otherTemperatures.length; i++) {
            double fac = Math.sqrt(otherTemperatures[i]/temperature);
            
            for (int j=0; j<atoms.getAtomCount(); j++) {
                IAtom jRealAtom = atoms.getAtom(j);
                IVectorMutable pos = pretendAtoms.getAtom(j).getPosition();
                pos.Ea1Tv1(1-fac, coordinateDefinition.getLatticePosition(jRealAtom));
                pos.PEa1Tv1(+fac, jRealAtom.getPosition());
            }

            double otherU = meterPotential.getDataAsScalar();
            double ai = (otherU-latticeEnergy)/otherTemperatures[i];
//            if (i==0) System.out.println(u+" "+otherU+" "+(u-latticeEnergy)/temperature+" "+(otherU-latticeEnergy)/otherTemperatures[i]);
            
            for (int j=0; j<numAlpha; j++) {
                if (temperature>otherTemperatures[i]) {
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
    public static void openFW(String filename) {
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
    public static void closeFW() {
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

    public CoordinateDefinition getCoordinateDefinition() {
        return coordinateDefinition;
    }

    public void setCoordinateDefinition(CoordinateDefinition newCoordinateDefinition) {
        this.coordinateDefinition = newCoordinateDefinition;

        // insert atoms into the box at their lattice sites.
        // we do this because want to find neighbors now (and then never again)
        IBox realBox = coordinateDefinition.getBox();
        pretendBox.setBoundary(realBox.getBoundary());
        pretendBox.setNMolecules(species, realBox.getNMolecules(species));
        IAtomList atoms = realBox.getLeafList();
        IAtomList pretendAtoms = pretendBox.getLeafList();
        for (int j=0; j<atoms.getAtomCount(); j++) {
            IAtom jRealAtom = atoms.getAtom(j);
            IVectorMutable pos = pretendAtoms.getAtom(j).getPosition();
            pos.E(coordinateDefinition.getLatticePosition(jRealAtom));
        }

        if (potentialMaster instanceof PotentialMasterList) {
            // find neighbors now.
            ((PotentialMasterList)potentialMaster).getNeighborManager(pretendBox).reset();
        }
    }

}
