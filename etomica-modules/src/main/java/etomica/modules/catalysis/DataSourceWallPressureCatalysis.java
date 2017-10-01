/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.catalysis;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataSource;
import etomica.data.IDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorHard.CollisionListener;
import etomica.modules.catalysis.InteractionTracker.CatalysisAgent;
import etomica.potential.P1HardBoundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.units.dimensions.Pressure;

public class DataSourceWallPressureCatalysis implements IDataSource, CollisionListener {
    public DataSourceWallPressureCatalysis(Space space, ISpecies speciesC, ISpecies speciesO, AtomLeafAgentManager interactionAgentManager) {
        this.space = space;
        this.speciesC = speciesC;
        this.speciesO = speciesO;
        this.interactionAgentManager = interactionAgentManager;
        data = new DataDoubleArray(3);
        dataInfo = new DataInfoDoubleArray("pressure", Pressure.DIMENSION, new int[]{3});
        tag = new DataTag();
        dataInfo.addTag(tag);
    }

    /**
     * Implementation of CollisionListener interface
     * Adds collision virial (from potential) to accumulator
     */
    public void collisionAction(IntegratorHard.Agent agent) {
        if (agent.collisionPotential instanceof P1HardBoundary) {
            Vector p = agent.atom.getPosition();
            if (p.getX(1) > 0) {
                IAtom atom = agent.atom;
                CatalysisAgent catalysisAgent = (CatalysisAgent)interactionAgentManager.getAgent(atom);
                if (atom.getType() == speciesC.getAtomType(0)) {
                    if (catalysisAgent.bondedAtom2 == null) {
                        virialSumCO += ((P1HardBoundary)agent.collisionPotential).lastWallVirial();
                    }
                    else {
                        virialSumCO2 += ((P1HardBoundary)agent.collisionPotential).lastWallVirial();
                    }
                }
                else if (atom.getType() == speciesO.getAtomType(0)) {
                    if (catalysisAgent.bondedAtom1.getType() == atom.getType()) {
                        virialSumO2 += ((P1HardBoundary)agent.collisionPotential).lastWallVirial();
                    }
                    else {
                        CatalysisAgent catalysisAgentC = (CatalysisAgent)interactionAgentManager.getAgent(catalysisAgent.bondedAtom1);
                        if (catalysisAgentC.bondedAtom2 == null) {
                            virialSumCO += ((P1HardBoundary)agent.collisionPotential).lastWallVirial();
                        }
                        else {
                            virialSumCO2 += ((P1HardBoundary)agent.collisionPotential).lastWallVirial();
                        }
                    }

                }
                else {
                    throw new RuntimeException("oops");
                }
            }
        }
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    public IData getData() {
        double currentTime = integratorHard.getCurrentTime();
        double[] x = data.getData();
        x[0] = virialSumCO / (currentTime - lastTime);
        x[1] = virialSumO2 / (currentTime - lastTime);
        x[2] = virialSumCO2 / (currentTime - lastTime);
        lastTime = currentTime;
        virialSumCO = virialSumO2 = virialSumCO2 = 0;
        return data;
    }

    /**
     * Registers meter as a collisionListener to the integrator, and sets up
     * a DataSourceTimer to keep track of elapsed time of integrator.
     */
    public void setIntegrator(IntegratorHard newIntegrator) {
        if(newIntegrator == integratorHard) return;
        if(integratorHard != null) {
            integratorHard.removeCollisionListener(this);
        }
        integratorHard = newIntegrator;
        if(newIntegrator != null) {
            integratorHard.addCollisionListener(this);
            lastTime = integratorHard.getCurrentTime();
        }
        virialSumCO = virialSumO2 = virialSumCO2 = 0;
    }

    public IntegratorHard getIntegrator() {
        return integratorHard;
    }

    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected final ISpecies speciesC, speciesO;
    protected final AtomLeafAgentManager interactionAgentManager;
    protected IntegratorHard integratorHard;
    protected double virialSumCO, virialSumO2, virialSumCO2;
    protected double lastTime;
    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
}
