/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.space.Vector;
import etomica.atom.iterator.IteratorDirective;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.nbr.list.PotentialMasterList;
import etomica.space.Space;
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

    public MeterSolidDACut(Space space, PotentialMaster potentialMaster, CoordinateDefinition coordinateDefinition, double[] cutoffs) {
        this.coordinteDefinition = coordinateDefinition;
        tag = new DataTag();
        this.potentialMaster = potentialMaster;
        iteratorDirective = new IteratorDirective();
        iteratorDirective.includeLrc = false;
        dim = space.D();
        box = coordinateDefinition.getBox();

        latticeEnergyDADv2 = latticeEnergy = new double[cutoffs.length];
        latticePressureDADv2 = latticePressure = new double[cutoffs.length];

        if (potentialMaster instanceof PotentialMasterList) {
            pc = new PotentialCalculationSolidSuperCut(space, coordinateDefinition, cutoffs);
        }
        else {
            pc = new PotentialCalculationSolidSuperCutLS(space, coordinateDefinition, cutoffs);
        }
        pc.zeroSum();
        potentialMaster.calculate(box, iteratorDirective, pc);
        double[] energy = pc.getEnergySum();
        double[] virial = pc.getVirialSum();
        
        if (false) System.out.print("Lattice energy: ");
        for (int i=0; i<cutoffs.length; i++) {
            latticeEnergy[i] = energy[i];
            if (false) System.out.print(" "+latticeEnergy[i]/box.getMoleculeList().getMoleculeCount());
            latticePressure[i] = -virial[i]/(box.getBoundary().volume()*dim);
        }
        if (false) System.out.print("\n");
        if (false) {
            System.out.print("Lattice pressure: ");
            for  (int i=0; i<cutoffs.length; i++) {
                System.out.print(" "+latticePressure[i]);
            }
            System.out.println();
        }

        int n = 6*cutoffs.length;
        dataInfo = new DataInfoDoubleArray("Stuff", Null.DIMENSION, new int[]{n});
        dataInfo.addTag(tag);
        data = new DataDoubleArray(n);

        if (potentialMaster instanceof PotentialMasterList) {
            pcDADv2 = new PotentialCalculationSolidSuperCut(space, coordinateDefinition, cutoffs);
        }
        else {
            pcDADv2 = new PotentialCalculationSolidSuperCutLS(space, coordinateDefinition, cutoffs);
        }
    }

    public void setPotentialMasterDADv2(PotentialMaster potentialMasterDADv2, double[] bpResDADv2) {
        this.potentialMasterDADv2 = potentialMasterDADv2;
        this.bpResDADv2 = bpResDADv2;
        latticeEnergyDADv2 = new double[latticeEnergy.length];
        latticePressureDADv2 = new double[latticeEnergy.length];

        pcDADv2.zeroSum();
        potentialMasterDADv2.calculate(box, iteratorDirective, pcDADv2);
        double[] energy = pcDADv2.getEnergySum();
        double[] virial = pcDADv2.getVirialSum();
        if (false) System.out.print("LJ Lattice energy: ");
        for (int i=0; i<energy.length; i++) {
            latticeEnergyDADv2[i] = energy[i];
            if (false) System.out.print(" "+latticeEnergyDADv2[i]/box.getMoleculeList().getMoleculeCount());
            latticePressureDADv2[i] = -virial[i]/(box.getBoundary().volume()*dim);
        }
        if (false) {
            System.out.println();
            System.out.print("LJ Lattice pressure: ");
            for  (int i=0; i<latticePressureDADv2.length; i++) {
                System.out.print(" "+latticePressureDADv2[i]);
            }
            System.out.println();
        }
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
        bpResDADv2 = bpRes;
    }
    
    /**
     * Computes total pressure in box by summing virial over all pairs, and adding
     * ideal-gas contribution.
     */
    public IData getData() {
        double V = box.getBoundary().volume();
        int N = box.getMoleculeList().getMoleculeCount();
        double rho = N/V;
    	pc.zeroSum();
        potentialMaster.calculate(box, iteratorDirective, pc);
        double[] p1 = pc.getPressure1();
        Vector[] p1XYZ = pc.getPressure1XYZ();
        double[] virial = pc.getVirialSum();
        double[] energy = pc.getEnergySum();
        double[] dadb = pc.getDADBSum();
        Vector[] dadbXYZ = pc.getDADBXYZ();

        double[] p1DADv2 = p1;
        double[] energyDADv2 = energy;
        double[] dadbDADv2 = dadb;

        if (potentialMasterDADv2 != null) {
            pcDADv2.zeroSum();
            potentialMasterDADv2.calculate(box, iteratorDirective, pcDADv2);
            p1DADv2 = pcDADv2.getPressure1();
            energyDADv2 = pcDADv2.getEnergySum();
            dadbDADv2 = pcDADv2.getDADBSum();
        }

        double[] x = data.getData();
        int j=0;
        for (int i=0; i<p1.length; i++) {
            double measuredP = temperature*rho - virial[i]/(dim*V);
            double buc = (0.5*dadb[i] + (energy[i] - latticeEnergy[i]))/temperature/N;
            double fac2 = (-1/V + bpRes[i])/(dim*N-dim);
            double Zc = (p1[i] + fac2*dadb[i] - latticePressure[i])/(rho*temperature);
            double PcZaniso = 3*((p1XYZ[i].getX(2) - 0.5*(p1XYZ[i].getX(0)+p1XYZ[i].getX(1))) + fac2*(dadbXYZ[i].getX(2) - 0.5*(dadbXYZ[i].getX(0)+dadbXYZ[i].getX(1))));
            x[j+0] = energy[i]/N;
            x[j+1] = measuredP;
            x[j+2] = buc;
            // P = Plat + Pres + x[5]
            x[j+3] = Zc;
            // Pc = x[5]
            // Zc = x[5] / (rho*T)
            // this is dbAc/dv2 at constant Y (for LJ)
            double fac3 = 0;
            if (potentialMasterDADv2 != null) {
                buc = (0.5*dadbDADv2[i] + (energyDADv2[i] - latticeEnergyDADv2[i]))/temperature/N;
                fac3 = (bpResDADv2[i])/(dim*N-dim);
                Zc = (p1DADv2[i] + fac2*dadbDADv2[i] + fac3*dadb[i] - latticePressureDADv2[i])/(rho*temperature);
            }

            x[j+4] = (4*buc-Zc)*rho*rho/2;
            
            x[j+5] = PcZaniso;
            j+=6;
        }

        return data;
    }

    protected final int dim;
    protected final DataTag tag;
    protected DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected final PotentialMaster potentialMaster;
    protected PotentialMaster potentialMasterDADv2;
    private IteratorDirective iteratorDirective;
    private final PotentialCalculationSolidSuperCut pc, pcDADv2;
    protected double temperature;
    protected double[] latticeEnergy, latticePressure;
    protected double[] latticeEnergyDADv2, latticePressureDADv2;
    protected final Box box;
    protected double[] bpRes, bpResDADv2;
    protected final CoordinateDefinition coordinteDefinition;
}
