/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.atom.iterator.IteratorDirective;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
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
 
public class MeterSolidDACut implements IEtomicaDataSource {

    public MeterSolidDACut(ISpace space, IPotentialMaster potentialMaster, CoordinateDefinition coordinateDefinition, double[] cutoffs) {
        this.coordinteDefinition = coordinateDefinition;
        tag = new DataTag();
        this.potentialMaster = potentialMaster;
        iteratorDirective = new IteratorDirective();
        iteratorDirective.includeLrc = false;
        dim = space.D();
        box = coordinateDefinition.getBox();

        latticeEnergy = new double[cutoffs.length];
        latticePressure = new double[cutoffs.length];
        pc = new PotentialCalculationSolidSuperCut(space, coordinateDefinition, cutoffs);
        pc.zeroSum();
        potentialMaster.calculate(box, iteratorDirective, pc);
        double[] energy = pc.getEnergySum();
        double[] virial = pc.getVirialSum();
        System.out.print("Lattice energy: ");
        for (int i=0; i<cutoffs.length; i++) {
            latticeEnergy[i] = energy[i];
            System.out.print(" "+latticeEnergy[i]/box.getMoleculeList().getMoleculeCount());
            latticePressure[i] = -virial[i]/(box.getBoundary().volume()*dim);
        }
        System.out.println();
        System.out.print("Lattice pressure: ");
        for  (int i=0; i<cutoffs.length; i++) {
            System.out.print(" "+latticePressure[i]);
        }
        System.out.println();

        int n = 5*cutoffs.length;
        dataInfo = new DataInfoDoubleArray("Stuff", Null.DIMENSION, new int[]{n});
        dataInfo.addTag(tag);
        data = new DataDoubleArray(n);
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
    
    public void setBPRes(double[] bpRes) {
        this.bpRes = bpRes;
    }
    
    /**
     * Computes total pressure in box by summing virial over all pairs, and adding
     * ideal-gas contribution.
     */
    public IData getData() {
    	pc.zeroSum();
        potentialMaster.calculate(box, iteratorDirective, pc);
        double V = box.getBoundary().volume();
        int N = box.getMoleculeList().getMoleculeCount();
        double rho = N/V;
        double[] x = data.getData();
        double[] p1 = pc.getPressure1();
        double[] virial = pc.getVirialSum();
        double[] energy = pc.getEnergySum();
        double[] dadb = pc.getDADBSum();
        
        int j=0;
        for (int i=0; i<p1.length; i++) {
            double measuredP = temperature*rho - virial[i]/(dim*V);
            x[j+0] = energy[i]/N;
            x[j+1] = measuredP;
            double buc = (0.5*dadb[i] + (energy[i] - latticeEnergy[i]))/temperature/N;
            x[j+2] = buc;
            double fac2 = (-1/V + bpRes[i])/(dim*N-dim);
            // P = Plat + Pres + x[5]
            double Zc = (p1[i] + fac2*dadb[i] - latticePressure[i])/(rho*temperature);
            x[j+3] = Zc;
            // Pc = x[5]
            // Zc = x[5] / (rho*T)
            // this is dbAc/dv2 at constant Y (for LJ)
            x[j+4] = (4*buc-Zc)*rho*rho/2;
            j+=5;
        }

        return data;
    }

    protected final int dim;
    protected final DataTag tag;
    protected DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected final IPotentialMaster potentialMaster;
    private IteratorDirective iteratorDirective;
    private final PotentialCalculationSolidSuperCut pc;
    protected double temperature;
    protected double[] latticeEnergy, latticePressure;
    protected final IBox box;
    protected double[] bpRes;
    protected final CoordinateDefinition coordinteDefinition;
}
