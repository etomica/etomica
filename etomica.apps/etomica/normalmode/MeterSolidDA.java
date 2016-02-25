/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.IVectorMutable;
import etomica.atom.iterator.IteratorDirective;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.potential.PotentialCalculationEnergySum;
import etomica.potential.PotentialCalculationVirialSum;
import etomica.space.ISpace;
import etomica.units.Null;

/**
 * Meter for evaluation of the soft-potential pressure in a box.
 * Requires that temperature be set in order to calculation ideal-gas
 * contribution to pressure; default is to use zero temperature, which
 * causes this contribution to be omitted.
 *
 * @author David Kofke
 */
 
public class MeterSolidDA implements IEtomicaDataSource {

    public MeterSolidDA(ISpace space, IPotentialMaster potentialMaster, CoordinateDefinition coordinateDefinition, boolean doD2) {
        this.coordinteDefinition = coordinateDefinition;
        tag = new DataTag();
    	this.potentialMaster = potentialMaster;
        iteratorDirective = new IteratorDirective();
        iteratorDirective.includeLrc = false;
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
        

        int n = 5;
        dataInfo = new DataInfoDoubleArray("Stuff", Null.DIMENSION, new int[]{n});
        dataInfo.addTag(tag);
        data = new DataDoubleArray(n);
        this.doD2 = doD2;
        dr = space.makeVector();
    }
    
    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public void setTemperature(double temperature) {
        this.temperature = temperature;
    }
    
    public void setPRes(double pRes) {
        this.pRes = pRes;
    }
    
    /**
     * Computes total pressure in box by summing virial over all pairs, and adding
     * ideal-gas contribution.
     */
    public IData getData() {
    	pc.zeroSum();
        potentialMaster.calculate(box, iteratorDirective, pc);
        double p1 = pc.getPressure1();
        double[] x = data.getData();
        double V = box.getBoundary().volume();
        double rho = box.getMoleculeList().getMoleculeCount()/V;
        double measuredP = temperature*rho - pc.getVirialSum()/(dim*V);
        int N = box.getMoleculeList().getMoleculeCount();
        x[0] = pc.getEnergySum()/N;
        x[1] = measuredP;
        double buc = (0.5*pc.getDADBSum() + (pc.getEnergySum() - latticeEnergy))/temperature/N;
        x[2] = buc;
        double vol = box.getBoundary().volume();
        double fac2 = (-1/vol + pRes/temperature)/(dim*N-dim);
        // P = Plat + Pres + x[5]
        double density = N / vol;
        double Zc = (p1 + fac2*pc.getDADBSum() - latticePressure)/(density*temperature);
        x[3] = Zc;
        // Pc = x[5]
        // Zc = x[5] / (rho*T)
        // this is dbAc/dv2 at constant Y (for LJ)
        x[4] = (4*buc-Zc)*density*density/2;

        return data;
    }

    protected final int dim;
    protected final DataTag tag;
    protected DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected final IPotentialMaster potentialMaster;
    private IteratorDirective iteratorDirective;
    private final PotentialCalculationSolidSuper pc;
    protected double temperature;
    protected double latticeEnergy, latticePressure;
    protected final IBox box;
    protected double pRes;
    protected final boolean doD2;
    protected final CoordinateDefinition coordinteDefinition;
    protected final IVectorMutable dr;
}
