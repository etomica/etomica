/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.rosmosis;

import etomica.atom.MoleculePositionGeometricCenter;
import etomica.box.Box;
import etomica.atom.IMolecule;
import etomica.atom.IMoleculeList;
import etomica.simulation.Simulation;
import etomica.api.ISpecies;
import etomica.space.Vector;
import etomica.atom.IMoleculePositionDefinition;
import etomica.atom.MoleculeAgentManager;
import etomica.atom.MoleculeAgentManager.MoleculeAgentSource;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorMD;
import etomica.space.Space;
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
public class MeterFlux implements IEtomicaDataSource, MoleculeAgentSource {

    public MeterFlux(Simulation sim, Space _space) {
        this.sim = sim;
    	this.space = _space;
        data = new DataDouble();
        dataInfo = new DataInfoDouble("flux", new CompoundDimension(new Dimension[]{
                Quantity.DIMENSION, Time.DIMENSION, Length.DIMENSION}, new double[]{1,-1,0}));
        tag = new DataTag();
        boundaries = new double[0];
        positionDefinition = new MoleculePositionGeometricCenter(space);
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
            agentManager = new MoleculeAgentManager(sim, box, this);
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
            agentManager = new MoleculeAgentManager(sim, box, this);
        }
    }
    
    public Box getBox() {
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
    
    public IData getData() {
        int crossings = 0;
        double boxLength = box.getBoundary().getBoxSize().getX(dim);
        for (int i=0; i<species.length; i++) {
            IMoleculeList molecules = box.getMoleculeList(species[i]);
            for (int j=0; j<molecules.getMoleculeCount(); j++) {
                IMolecule atom = molecules.getMolecule(j);
                Vector oldPosition = ((Vector)agentManager.getAgent(atom));
                double oldX = oldPosition.getX(dim);
                Vector newPosition = positionDefinition.position(atom);
                double newX = newPosition.getX(dim);
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
        double newTime;
        if (integrator instanceof IntegratorMD) {
            newTime = ((IntegratorMD)integrator).getCurrentTime();
            if (newTime < oldTime) {
                // reinitialize sets time back to 0
                oldTime = newTime - ((IntegratorMD)integrator).getTimeStep();
            }
        }
        else {
            newTime = integrator.getStepCount();
            if (newTime < oldTime) {
                oldTime = newTime - 1;
            }
        }
        data.x /= newTime - oldTime;
        oldTime = newTime;
        for (int i=0; i<space.D(); i++) {
            if (i == dim) continue;
            data.x /= box.getBoundary().getBoxSize().getX(i);
        }
        return data;
    }

    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    public Class getMoleculeAgentClass() {
        return Vector.class;
    }

    public Object makeAgent(IMolecule a) {
        ISpecies thisSpecies = a.getType();
        for (int i=0; i<species.length; i++) {
            if (species[i] == thisSpecies) {
                Vector vec = space.makeVector();
                vec.E(positionDefinition.position(a));
                return vec;
            }
        }
        // not one of the species we care about
        return null;
    }

    public void releaseAgent(Object agent, IMolecule atom) {
        /* do nothing */
    }

    protected final Simulation sim;
    protected final DataDouble data;
    protected DataInfoDouble dataInfo;
    protected final DataTag tag;
    protected ISpecies[] species;
    protected Box box;
    protected double[] boundaries;
    protected int[] boundaryCoefficients;
    protected int dim;
    protected MoleculeAgentManager agentManager;
    protected IntegratorBox integrator;
    protected double oldTime;
    protected long oldStep;
    private final Space space;
    protected IMoleculePositionDefinition positionDefinition;
}
