package etomica.modules.rosmosis;

import etomica.api.IBox;
import etomica.api.IVector;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.box.Box;
import etomica.data.Data;
import etomica.data.DataSource;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntegratorNonintervalEvent;
import etomica.integrator.IntegratorNonintervalListener;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.units.CompoundDimension;
import etomica.units.Dimension;
import etomica.units.Length;
import etomica.units.Quantity;
import etomica.units.Time;

/**
 * Meter to measure flux across a boundary or boundaries.  If an atom is on one
 * side of the boundary at one time and on the other side of the boundary the
 * next time, then the meter counts that as a crossing.  Each boundary can have
 * a coefficient, such that flow into one region through boundaries on either
 * side can be considered "positive" flux.
 * 
 * If an IntegratorMD is used, flux will be given in terms of crossings per
 * area per time.  Otherwise, flux will be given in terms of crossings per
 * area per step.
 *
 * @author Andrew Schultz
 */
public class MeterFlux implements DataSource, AgentSource, IntegratorNonintervalListener {

    public MeterFlux(Space _space) {
    	this.space = _space;
        data = new DataDouble();
        dataInfo = new DataInfoDouble("flux", new CompoundDimension(new Dimension[]{
                Quantity.DIMENSION, Time.DIMENSION, Length.DIMENSION}, new double[]{1,-1,0}));
        tag = new DataTag();
        boundaries = new double[0];
    }
    
    public void setBoundaries(int newDim, double[] newBoundaries, int[] newBoundaryCoefficients) {
        dim = newDim;
        boundaries = newBoundaries;
        boundaryCoefficients = newBoundaryCoefficients;
    }
    
    public int getDim() {
        return dim;
    }
    
    public double[] getBoundaries() {
        return boundaries;
    }
    
    public void setSpecies(ISpecies[] newSpecies) {
        species = newSpecies;
        if (box != null) {
            agentManager = new AtomAgentManager(this, box);
        }
    }
    
    public ISpecies[] getSpecies() {
        return species;
    }
    
    public void setBox(Box newBox) {
        box = newBox;
        if (integrator != null) {
            if (integrator instanceof IntegratorMD) {
                dataInfo = new DataInfoDouble("flux", new CompoundDimension(new Dimension[]{
                        Quantity.DIMENSION, Time.DIMENSION, Length.DIMENSION}, new double[]{1,-1,1-space.D()}));
            }
            else {
                // Quantity(crossings) / Quantity(steps)
                dataInfo = new DataInfoDouble("flux", new CompoundDimension(new Dimension[]{
                        Length.DIMENSION}, new double[]{1-space.D()}));
            }
        }
        if (species != null) {
            agentManager = new AtomAgentManager(this, box);
        }
    }
    
    public IBox getBox() {
        return box;
    }
    
    public void setIntegrator(IntegratorBox newIntegrator) {
        integrator = newIntegrator;
        if (integrator instanceof IntegratorMD) {
            oldTime = ((IntegratorMD)integrator).getCurrentTime();
            if (box != null) {
                dataInfo = new DataInfoDouble("flux", new CompoundDimension(new Dimension[]{
                        Quantity.DIMENSION, Time.DIMENSION, Length.DIMENSION}, new double[]{1,-1,1-space.D()}));
            }
        }
        else {
            oldStep = integrator.getStepCount();
            if (box != null) {
                // Quantity(crossings) / Quantity(steps)
                dataInfo = new DataInfoDouble("flux", new CompoundDimension(new Dimension[]{
                        Length.DIMENSION}, new double[]{1-space.D()}));
            }
        }
    }
    
    public Data getData() {
        int crossings = 0;
        double boxLength = box.getBoundary().getDimensions().x(dim);
        for (int i=0; i<species.length; i++) {
            AtomSet molecules = box.getMoleculeList(species[i]);
            for (int j=0; j<molecules.getAtomCount(); j++) {
                IAtom atom = molecules.getAtom(j);
                IVector oldPosition = ((IVector)agentManager.getAgent(atom));
                double oldX = oldPosition.x(dim);
                IVector newPosition = atom.getType().getPositionDefinition().position(atom);
                double newX = newPosition.x(dim);
                for (int k=0; k<boundaries.length; k++) {
                    double newDelta = newX - boundaries[k];
                    if (Math.abs(newDelta)  > 0.25*boxLength) continue;
                    double oldDelta = oldX - boundaries[k];
                    if (Math.abs(oldDelta)  > 0.25*boxLength) continue;
                    if (newDelta * oldDelta < 0) {
                        crossings += boundaryCoefficients[k] * Math.abs(newDelta)/newDelta;
                        break;
                    }
                }
                oldPosition.E(newPosition);
            }
        }
        data.x = crossings;
        if (integrator instanceof IntegratorMD) {
            data.x /= ((IntegratorMD)integrator).getCurrentTime() - oldTime;
            oldTime = ((IntegratorMD)integrator).getCurrentTime();
        }
        else {
            data.x /= integrator.getStepCount() - oldStep;
            oldTime = ((IntegratorMD)integrator).getCurrentTime();
        }
        for (int i=0; i<space.D(); i++) {
            if (i == dim) continue;
            data.x /= box.getBoundary().getDimensions().x(i);
        }
        return data;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    public Class getAgentClass() {
        return IVector.class;
    }

    public Object makeAgent(IAtom a) {
        if (a instanceof IAtomPositioned) {
            // oh, the irony
            return null;
        }
        ISpecies thisSpecies = a.getType().getSpecies();
        for (int i=0; i<species.length; i++) {
            if (species[i] == thisSpecies) {
                IVector vec = space.makeVector();
                vec.E(a.getType().getPositionDefinition().position(a));
                return vec;
            }
        }
        // not one of the species we care about
        return null;
    }

    public void releaseAgent(Object agent, IAtom atom) {
        /* do nothing */
    }

    public void nonintervalAction(IntegratorNonintervalEvent evt) {
        if (evt.type() == IntegratorNonintervalEvent.RESET &&
                agentManager != null) {
            agentManager = new AtomAgentManager(this, box);
        }
        oldTime = 0;
    }

    protected final DataDouble data;
    protected DataInfoDouble dataInfo;
    protected final DataTag tag;
    protected ISpecies[] species;
    protected Box box;
    protected double[] boundaries;
    protected int[] boundaryCoefficients;
    protected int dim;
    protected AtomAgentManager agentManager;
    protected IntegratorBox integrator;
    protected double oldTime;
    protected long oldStep;
    private final Space space;
}
