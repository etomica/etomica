package etomica;
import etomica.units.Dimension;

/**
 * Takes measurments to average pressure from momentum flux at Boundary of Space
 * Intended for use with BoundaryHard
 */
 
public class MeterPressureBoundary extends Meter
{
    private double momentumSum = 0.0;
    private double timeSum = 0.0;

    public MeterPressureBoundary() {
        this(Simulation.instance);
    }
    public MeterPressureBoundary(Simulation sim)
    {
        super(sim);
        setLabel("Pressure");
    }
    
    /**
     * Declaration that this meter does not use the boundary object of phase when making its measurements
     */
    protected final boolean usesPhaseBoundary() {return false;}
    /**
     * Declaration that this meter does not use the iteratorFactory of phase when making its measurements
     */
    protected final boolean usesPhaseIteratorFactory() {return false;}
        
    public Dimension getDimension() {return Dimension.PRESSURE;}
    
    public void intervalAction(Integrator.IntervalEvent evt) {
        IntegratorMD integrator = (IntegratorMD)evt.getSource();
        timeSum += integrator.timeStep * integrator.interval;
        updateSums();}

    public double currentValue()
    {
        Space2D.BoundaryHard bnd = (Space2D.BoundaryHard)phase.boundary();
        double flux = bnd.pAccumulator/timeSum;   //divide by time interval
        flux /= (2*(bnd.dimensions.x+bnd.dimensions.y)); //divide by area
        timeSum = 0.0;          //zeroing should be moved to intervalAction?
        bnd.pAccumulator = 0.0;
        return flux;
    }
}