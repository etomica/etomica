/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import java.io.FileWriter;
import java.io.IOException;
import java.lang.Math;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.potential.*;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Null;


/**
 * Meter for evaluation of the soft-potential pressure in a box.
 * Requires that temperature be set in order to calculation ideal-gas
 * contribution to pressure; default is to use zero temperature, which
 * causes this contribution to be omitted.
 *
 * @author David Kofke
 */

public class MeterPressureHMA implements IDataSource {
    protected final int dim;
    protected final DataTag tag;
    protected DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected final PotentialMaster potentialMaster;
    private IteratorDirective iteratorDirective;
    private final PotentialCalculationSolidSuper pc;
    private int N;
    private double latticeConstant;
    protected double temperature, truncationRadius, top, bottom;
    protected double latticeEnergy, latticePressure;
    protected double[] GruneisenParameterList, Gij;
    protected final PotentialCalculationForceSum pcForce;
    protected final Box box;
    protected double betaHarmonicPressure, ratio;
    protected final boolean doD2;
    protected final CoordinateDefinition coordinateDefinition;
    protected final AtomLeafAgentManager<Vector> forceManager;
    protected Vector dr;
    private FileWriter fileWriter;

    public MeterPressureHMA(Space space, PotentialMaster potentialMaster, CoordinateDefinition coordinateDefinition, boolean doD2) {
        this.coordinateDefinition = coordinateDefinition;
        this.potentialMaster = potentialMaster;
        this.doD2 = doD2;
        iteratorDirective = new IteratorDirective();
        iteratorDirective.includeLrc = false;
        tag = new DataTag();

        dim = space.D();
        box = coordinateDefinition.getBox();
        N = box.getMoleculeList().size();
        latticeConstant = box.getBoundary().volume() / N;

        PotentialCalculationEnergySum pcEnergy = new PotentialCalculationEnergySum();
        pcEnergy.zeroSum();
        potentialMaster.calculate(box, iteratorDirective, pcEnergy);
        latticeEnergy = pcEnergy.getSum();

        PotentialCalculationVirialSum pcVirial = new PotentialCalculationVirialSum();
        pcVirial.zeroSum();
        potentialMaster.calculate(box, iteratorDirective, pcVirial);
        latticePressure = -pcVirial.getSum() / (box.getBoundary().volume() * dim);

        pcForce = new PotentialCalculationForceSum();
        forceManager = new AtomLeafAgentManager<>(a -> space.makeVector(), box);
        pcForce.setAgentManager(forceManager);

        pc = new PotentialCalculationSolidSuper(space, coordinateDefinition);
        pc.setDoSecondDerivative(doD2);

        int n = doD2 ? 7 : 5;       // If doD2 is true, n = 7, else n = 5.
        dataInfo = new DataInfoDoubleArray("Stuff", Null.DIMENSION, new int[]{n});
        dataInfo.addTag(tag);
        data = new DataDoubleArray(n);
        dr = space.makeVector();

    }

    public void setTemperature(double temperature) {
        this.temperature = temperature;
    }

    // TODO: Can't we use the information from the simulation?
    public void setTruncationRadius(double truncationRadius) {
        this.truncationRadius = truncationRadius;
    }

    public void setPRes() {
        double[] GruneisenParameterList = calculateGruneisenParameters();
        betaHarmonicPressure = 0.0;
        for(int i=1; i < GruneisenParameterList.length; i++)
            betaHarmonicPressure += GruneisenParameterList[i];
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    public void calculateGij() {
        double PI = Math.PI;
        Gij = new double[N];

        GruneisenParameterList = calculateGruneisenParameters();
        if (N % 2 == 1) {
            for (int index=0; index <= N - 1; index++) {
                for (int k=1; k <= (N-1)/2; k++) {
                    Gij[index] += GruneisenParameterList[k] * Math.cos((2 * PI * k)/ N * index);
                }
                Gij[index] = 2 * Gij[index] / N;
            }
        }
        else {
            for(int index=0; index <= N - 1; index++) {
                for (int k=1; k <= N/2 - 1; k++) {
                    Gij[index] += GruneisenParameterList[k] * Math.cos((2 * PI * k)/ N * index);
                }
                Gij[index] = (Math.pow(-1, index) * GruneisenParameterList[N/2] + 2 * Gij[index]) / N;
            }
        }
    }

    private double[] calculateGruneisenParameters() {
        double[] GruneisenParameterList = new double [dim * N];
        double PI = Math.PI;
        int numNeighbors = calculateNearestNeighbors();

        for (int k=1; k <= dim * (N - 1); k++) {
            double numerator = 0.0;
            double denominator = 0.0;

            for (int m=1; m <= numNeighbors; m++) {
                double reciprocal = 1.0 / (m * latticeConstant);
                numerator = numerator + m * (8736 * Math.pow(reciprocal, 15) - 1344 * Math.pow(reciprocal, 9)) * Math.pow(Math.sin(PI / N * k * m), 2) / N;
                denominator = denominator + (-624 * Math.pow(reciprocal, 14) +  168 * Math.pow(reciprocal, 8)) * Math.pow(Math.sin(PI / N * k * m), 2);
            }
            GruneisenParameterList[k] = -0.5 * numerator / denominator;
        }
        return GruneisenParameterList;
    }

    private int calculateNearestNeighbors() {
        // TODO: Should an exception be added for the case where we obtain '0'?
        return (int) (truncationRadius / latticeConstant);
    }

    public IData getData() {
        pc.zeroSum();
        potentialMaster.calculate(box, iteratorDirective, pc);

        double[] x = data.getData();
        double V = box.getBoundary().volume();
        double rho = N / V;
        double pressureConventional = temperature * rho - pc.getVirialSum() / (dim * V);
        double uSum = pc.getEnergySum();

        double pressureQuasiHarmonic = betaHarmonicPressure * temperature;
        x[0] = uSum / N;
        x[1] = (pressureConventional  - latticePressure - pressureQuasiHarmonic) / (temperature);
        double buc = (0.5 * pc.getDADBSum() + (uSum - latticeEnergy)) / temperature / N;
        x[2] = buc;
        x[3] = pressureQuasiHarmonic + latticePressure;
        double fV = (betaHarmonicPressure + 1.0 / V - rho) / (dim * (N - 1));
        double Zc = (-pc.getVirialSum() / (dim * V) + fV * pc.getDADBSum() - latticePressure) / (rho * temperature);
        double pressureHMAFromNormalMode = pressureQuasiHarmonic + temperature / V  - pc.getVirialSum() / (dim * V) - pc.getDADBSum() / V + calculateNormalModeSum();
        x[4] = (pressureHMAFromNormalMode - latticePressure - pressureQuasiHarmonic) / (temperature);
        double pressureHMAFromRealSpace = pressureQuasiHarmonic + temperature / V - pc.getVirialSum() / (dim * V) + fV * pc.getDADBSum();
        x[5] = (pressureHMAFromRealSpace - latticePressure - pressureQuasiHarmonic) / (temperature);
        x[6] = (4 * buc - Zc) * rho * rho / 2;

        return data;
    }

    private double calculateNormalModeSum() {
        IAtomList atoms = box.getLeafList();
        pcForce.reset();
        potentialMaster.calculate(box, iteratorDirective, pcForce);

        double normalModeSum = 0.0;
        for(int i=0; i < N; i++) {
            Vector forceActingOnAtom = forceManager.getAgent(atoms.get(i));
            for(int j=0; j < N; j++) {
                Vector xj = coordinateDefinition.getLatticePosition(atoms.get(j));
                Vector dr = box.getSpace().makeVector();
                dr.Ev1Mv2(atoms.get(j).getPosition(), xj);
                normalModeSum += forceActingOnAtom.dot(dr) * Gij[Math.abs(i - j)];
            }
        }
        return normalModeSum;
    }
}