package etomica.data.meter;
import etomica.EtomicaInfo;
import etomica.action.BoxInflate;
import etomica.action.BoxInflateDeformable;
import etomica.atom.iterator.IteratorDirective;
import etomica.data.Data;
import etomica.data.DataSource;
import etomica.data.DataSourceUniform;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.integrator.IntegratorBox;
import etomica.box.Box;
import etomica.potential.PotentialCalculationEnergySum;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.units.Pressure;
import etomica.units.Volume;

/**
 * Evaluates the pressure by examining the change in energy accompanying
 * small changes in volume.
 */
public class MeterPressureByVolumeChange implements DataSource, java.io.Serializable {
    
    public MeterPressureByVolumeChange(Space space) {
        this(space, makeDefaultDimensions(space.D()));
    }
    
    private static final boolean[] makeDefaultDimensions(int D) {
        boolean[] dim = new boolean[D];
        for (int i=0; i<D; i++) {
            dim[i] = true;
        }
        return dim;
    }
    
    public MeterPressureByVolumeChange(Space space, boolean[] dimensions) {
        this.space = space;
        tag = new DataTag();
        setX(-0.001, 0.001, 10);
        inflateDimensions = new boolean[space.D()];
        setInflateDimensions(dimensions);
        iteratorDirective = new IteratorDirective();
        inflater = new BoxInflate(space);
    }
    
    public MeterPressureByVolumeChange(Space space, BoxInflateDeformable pid){
        this(space, makeDefaultDimensions(space.D()), pid);
    }
    
    public MeterPressureByVolumeChange(Space space, boolean[] dimensions, BoxInflateDeformable pid){
        this.space = space;
        tag = new DataTag();
        setX(-0.001, 0.001, 10);
        inflateDimensions = new boolean[space.D()];
        setInflateDimensions(dimensions);
        iteratorDirective = new IteratorDirective();
        inflater = pid;
        
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Pressure measured by trial volume perturbations");
        return info;
    }
    
    /**
     * Sets the integrator associated with this instance.  The pressure is 
     * calculated for the box the integrator acts on and the integrator's 
     * temperature is used to calculate the pressure.
     */
    public void setIntegrator(IntegratorBox newIntegrator) {
        integrator = newIntegrator;
    }
    
    /**
     * Returns the integrator associated with this instance.  The pressure is 
     * calculated for the box the integrator acts on and the integrator's 
     * temperature is used to calculate the pressure.
     */
    public IntegratorBox getIntegrator() {
        return integrator;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }

    /**
     * For anisotropic volume change, indicates dimension in which volume is perturbed.
     */
    public final void setInflateDimensions(boolean[] directions) {
        if(directions.length != space.D()) {
            throw new IllegalArgumentException();
        }
        nDimension = 0;
        for(int i=0; i<directions.length; i++) {
            inflateDimensions[i] = directions[i];
            if(inflateDimensions[i]) nDimension++;
        }
        updateScale();
    }
    /**
     * Accessor method for setInflateDimension.
     */
    public boolean[] getInflateDimensions() {return inflateDimensions;}
    
    public void setX(double min, double max, int n) {
        xDataSource = new DataSourceUniform("x", Volume.dimension(space.D()), n, min, max);
        data = new DataDoubleArray(n);
        dataInfo = new DataInfoDoubleArray("Pressure by Volume Change", Pressure.dimension(space.D()), new int[]{n});
        dataInfo.addTag(tag);
        dataArray = data.getData();
        updateScale();
    }
    
    protected void updateScale() {
        if (inflateDimensions == null || xDataSource == null) {
            return;
        }
        
        //x is scaling in volume if isotropic, but is linear scaling if not isotropic
        double dx = (xDataSource.getXMax()-xDataSource.getXMin())/xDataSource.getNValues();
        final double[] x = ((DataDoubleArray)xDataSource.getData()).getData();
        for(int i=0; i<x.length; i++) { //disallow x = 0
            if(x[i] == 0.0) x[i] = 0.1*dx;
        }
        scale = new IVector[x.length];
        
        double mult = 1.0/nDimension;
        for(int i=0; i<x.length; i++) {
            scale[i] = space.makeVector();
            scale[i].E(Math.exp(mult*x[i]));
            for(int j=0; j<space.D(); j++) {
                if(!inflateDimensions[j]) scale[i].setX(j,1.0);
            }
//            System.out.println("scale " + i +" = " + scale[i]);
        }
    }
    
    public DataSource getXDataSource() {
        return xDataSource;
    }
    
    public Data getData() {
        if (integrator == null) throw new IllegalStateException("must call setIntegrator before using meter");
        Box box = integrator.getBox();
        inflater.setBox(box);
        energy.zeroSum();
        integrator.getPotential().calculate(box, iteratorDirective, energy);
        double uOld = energy.getSum();
        final double[] x = ((DataDoubleArray)xDataSource.getData()).getData();

        for(int i=0; i<dataArray.length; i++) {
            inflater.setVectorScale(scale[i]);
            inflater.actionPerformed();
            energy.zeroSum();
            integrator.getPotential().calculate(box, iteratorDirective, energy);
            double uNew = energy.getSum();
            dataArray[i] = Math.exp(-(uNew-uOld)/integrator.getTemperature()
                              + box.moleculeCount()*x[i]);
            inflater.undo();
        }
        

        return data;
    }

    private static final long serialVersionUID = 1L;
    private DataDoubleArray data;
    private IDataInfo dataInfo;
    private final DataTag tag;
    private double[] dataArray;
    private final BoxInflate inflater;
    private IVector[] scale;
    private final boolean[] inflateDimensions;
    private final IteratorDirective iteratorDirective;
    private final PotentialCalculationEnergySum energy = new PotentialCalculationEnergySum();
    private int nDimension;
    private final Space space;
    private DataSourceUniform xDataSource;
    private IntegratorBox integrator;
}
