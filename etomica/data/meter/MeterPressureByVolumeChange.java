package etomica.data.meter;
import etomica.EtomicaInfo;
import etomica.Phase;
import etomica.action.PhaseInflate;
import etomica.atom.iterator.IteratorDirective;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataSource;
import etomica.data.DataSourceUniform;
import etomica.data.types.DataDoubleArray;
import etomica.potential.PotentialCalculationEnergySum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.Dimension;

/**
 * Evaluates the pressure by examining the change in energy accompanying
 * small changes in volume.
 */
public class MeterPressureByVolumeChange implements Meter, java.io.Serializable {
    
    public MeterPressureByVolumeChange(PotentialMaster potentialMaster, boolean[] dimensions) {
        data = new DataDoubleArray("Pressure by Volume Change",Dimension.pressure(potentialMaster.getSpace().D()),10);
        dataArray = data.getData();
        spaceD = dimensions.length;
        potential = potentialMaster;
        setX(-0.001, 0.001, data.getLength());
        inflateDimensions = new boolean[spaceD];
        setInflateDimensions(dimensions);
        iteratorDirective = new IteratorDirective();
        inflater = new PhaseInflate(potentialMaster.getSpace());
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Pressure measured by trial volume perturbations");
        return info;
    }
    
    public DataInfo getDataInfo() {
        return data.getDataInfo();
    }

    /**
     * For anisotropic volume change, indicates dimension in which volume is perturbed.
     */
    public final void setInflateDimensions(boolean[] directions) {
        if(directions.length != spaceD) {
            throw new IllegalArgumentException();
        }
        nDimension = 0;
        for(int i=0; i<directions.length; i++) {
            inflateDimensions[i] = directions[i];
            if(inflateDimensions[i]) nDimension++;
            
        }
    }
    /**
     * Accessor method for setInflateDimension.
     */
    public boolean[] getInflateDimensions() {return inflateDimensions;}
    
    public void setX(double min, double max, int n) {
        xDataSource = new DataSourceUniform("x",Dimension.volume(potential.getSpace().D()),n, min, max);
        //x is scaling in volume if isotropic, but is linear scaling if not isotropic
        double dx = (max-min)/n;
        final double[] x = ((DataDoubleArray)xDataSource.getData()).getData();
        for(int i=0; i<x.length; i++) { //disallow x = 0
            if(x[i] == 0.0) x[i] = 0.1*dx;
        }
        scale = new Vector[x.length];
        
        double mult = 1.0/nDimension;
        for(int i=0; i<x.length; i++) {
            scale[i] = Space.makeVector(spaceD);
            scale[i].E(Math.exp(mult*x[i]));
            for(int j=0; j<spaceD; j++) {
                if(!inflateDimensions[j]) scale[i].setX(j,1.0);
            }
        }
    }
    
    /**
     * @return Returns the temperature.
     */
    public double getTemperature() {
        return temperature;
    }
    
    /**
     * @param temperature The temperature to set.
     */
    public void setTemperature(double temperature) {
        this.temperature = temperature;
    }
    
    public DataSource getXDataSource() {
        return xDataSource;
    }
    
    public Data getData() {
        if (phase == null) throw new IllegalStateException("must call setPhase before using meter");
        inflater.setPhase(phase);
        energy.zeroSum();
        potential.calculate(phase, iteratorDirective, energy);
        double uOld = energy.getSum();
        final double[] x = ((DataDoubleArray)xDataSource.getData()).getData();
        for(int i=0; i<dataArray.length; i++) {
            inflater.setVectorScale(scale[i]);
            inflater.actionPerformed();
            energy.zeroSum();
            potential.calculate(phase, iteratorDirective, energy);
            double uNew = energy.getSum();
            dataArray[i] = Math.exp(-(uNew-uOld)/temperature
                              + phase.moleculeCount()*x[i]);

            //TODO shouldn't this be done outside the loop?
            inflater.undo();
            //System.out.println( "  uNew " + uNew +" uOld " +uOld +" x " + x[i] +" scale" + scale[i]+ " y " +y[i] );
        }
        return data;
    }

    /**
     * @return Returns the phase.
     */
    public Phase getPhase() {
        return phase;
    }
    /**
     * @param phase The phase to set.
     */
    public void setPhase(Phase phase) {
        this.phase = phase;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }
    private String name;
    private Phase phase;
    private final DataDoubleArray data;
    private double[] dataArray;
    private final PhaseInflate inflater;
    Vector[] scale;
    boolean[] inflateDimensions;
    private IteratorDirective iteratorDirective;
    private final PotentialCalculationEnergySum energy = new PotentialCalculationEnergySum();
    private final PotentialMaster potential;
    private int nDimension;
    private int spaceD;
    private DataSourceUniform xDataSource;
    private double temperature = Double.NaN;
    
}
