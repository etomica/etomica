/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationEnergySum;
import etomica.potential.PotentialCalculationVirialSum;
import etomica.potential.PotentialMaster;
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

public class MeterSolidDA implements IDataSource {

    public MeterSolidDA(Space space, PotentialMaster potentialMaster, CoordinateDefinition coordinateDefinition, boolean doD2) {
        this.coordinteDefinition = coordinateDefinition;
        tag = new DataTag();
    	this.potentialMaster = potentialMaster;     // Passes in the potentialMaster.
        iteratorDirective = new IteratorDirective();
        iteratorDirective.includeLrc = false;       // Not including Long-range contributions.
        pc = new PotentialCalculationSolidSuper(space, coordinateDefinition);
        pc.setDoSecondDerivative(doD2);
        dim = space.D();
        box = coordinateDefinition.getBox();
        PotentialCalculationEnergySum pcEnergy = new PotentialCalculationEnergySum();
        pcEnergy.zeroSum();
        potentialMaster.calculate(box, iteratorDirective, pcEnergy);
        latticeEnergy = pcEnergy.getSum();

        PotentialCalculationVirialSum pcVirial = new PotentialCalculationVirialSum();
        pcVirial.zeroSum();
        potentialMaster.calculate(box, iteratorDirective, pcVirial);
        latticePressure = -pcVirial.getSum()/(box.getBoundary().volume()*dim);


        int n = doD2 ? 7 : 5;
        dataInfo = new DataInfoDoubleArray("Stuff", Null.DIMENSION, new int[]{n});
        dataInfo.addTag(tag);
        data = new DataDoubleArray(n);
        this.doD2 = doD2;
        dr = space.makeVector();
    }
    
    public IDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public void setTemperature(double temperature) {
        this.temperature = temperature;
        bpHarm = pRes / temperature;
    }
    
    public void setPRes(double pRes) {
        this.pRes = pRes;
        bpHarm = pRes / temperature;
    }
    
    /**
     * Computes total pressure in box by summing virial over all pairs, and adding
     * ideal-gas contribution.
     */
    public IData getData() {
    	pc.zeroSum();
        potentialMaster.calculate(box, iteratorDirective, pc);
        double[] x = data.getData();
        double V = box.getBoundary().volume();
        double rho = box.getMoleculeList().size()/V;
        double measuredP = temperature*rho - pc.getVirialSum()/(dim*V);
        int N = box.getMoleculeList().size();
        double uSum = pc.getEnergySum();
        x[0] = uSum / N;
        x[1] = measuredP;
        double buc = (0.5 * pc.getDADBSum() + (uSum - latticeEnergy)) / temperature / N;
        x[2] = buc;
        double vol = box.getBoundary().volume();
        // P = Plat + Pres + x[5]
        double density = N / vol;
        double fac2 = (-1 / V + pRes) / (dim * N - dim);
        double fV = (bpHarm - N / V) / (dim * (N - 1));
//        System.out.println(p1+" "+fac2+" "+pc.getDADBSum()+" "+latticePressure);
        double Zc = (-pc.getVirialSum() / (dim * V) + fV * pc.getDADBSum() - latticePressure) / (density * temperature);

        x[3] = Zc;
        // Pc = x[5]
        // Zc = x[5] / (rho*T)
        // this is dbAc/dv2 at constant Y (for LJ)


        x[4] = (4*buc-Zc)*density*density/2;
        //x[4] = (N-1)*1.5/N + (0.5 * pc.getDADBSum() + uSum) / temperature / N;

        if (doD2) {
            double y = latticeEnergy + dr.getD() * (N - 1) * temperature / 2;
            x[5] = (uSum - y) * (uSum - y);
//            System.out.println(uSum+" "+y+" "+x[5]);

            x[6] = -0.25 * (pc.getDADBSum() + pc.getD2Sum()) / temperature + (x[2] * x[2]) * N * N;
            //x[7] = -0.25 * (pc.getDADBSum() + pc.getD2Sum()) / temperature + x[4]*x[4]*N*N;
//            x[7] = pc.getD2Sum();
        }

        return data;
    }

    protected final int dim;
    protected final DataTag tag;
    protected DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected final PotentialMaster potentialMaster;
    private IteratorDirective iteratorDirective;
    private final PotentialCalculationSolidSuper pc;
    protected double temperature;
    protected double latticeEnergy, latticePressure;
    protected final Box box;
    protected double pRes, bpHarm;
    protected final boolean doD2;
    protected final CoordinateDefinition coordinteDefinition;
    protected final Vector dr;
}
