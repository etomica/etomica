/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode.nptdemo.fluid;

import etomica.action.BoxInflate;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.ISimulation;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.iterator.IteratorDirective;
import etomica.box.Box;
import etomica.data.DataSourceUniform;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.integrator.IntegratorBox;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.PotentialCalculationEnergySum;
import etomica.space.ISpace;
import etomica.units.Null;
import etomica.units.Pressure;
import etomica.units.Volume;
import etomica.util.Arrays;

/**
 * Evaluates the pressure by examining the change in energy accompanying
 * small changes in volume.
 *
 * Pressure is approximately log(<x>)/deltaV
 * where x is one of the values returned by getData (<x> is that value averaged)
 * and deltaV is the change in volume associated with that x value.  deltaV = s * V
 * where s is the scaling factor (available via getScalingDataSource)
 */
public class MeterPressureByVolumeChangeStuff implements IEtomicaDataSource, java.io.Serializable {
    
    public MeterPressureByVolumeChangeStuff(ISimulation sim, ISpace space, IStuff stuff) {
        this.sim = sim;
        this.space = space;
        tag = new DataTag();
        setX(-0.001, 0.001, 10);
        iteratorDirective = new IteratorDirective();
        inflater = new BoxInflate(space);
        myBox = new Box(space);
        sim.addBox(myBox);
        this.stuff = stuff;
    }
    
    /**
     * Sets the integrator associated with this instance.  The pressure is 
     * calculated for the box the integrator acts on and the integrator's 
     * temperature is used to calculate the pressure.
     */
    public void setIntegrator(IntegratorBox newIntegrator) {
        integrator = newIntegrator;
        myBox.setBoundary(integrator.getBox().getBoundary());
    }
    
    /**
     * Returns the integrator associated with this instance.  The pressure is 
     * calculated for the box the integrator acts on and the integrator's 
     * temperature is used to calculate the pressure.
     */
    public IntegratorBox getIntegrator() {
        return integrator;
    }

    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }

    public void setX(double min, double max, int n) {
        if (n < 1) throw new IllegalArgumentException("n must be greater than 1");
        xDataSource = new DataSourceUniform("x", Volume.dimension(space.D()), n, min, max);
        vDataSource = new DataSourceExp(xDataSource);
        vDataSource.reset();
        dataInfo = new DataInfoDoubleArray("Pressure by Volume Change", Pressure.dimension(space.D()), new int[]{vDataSource.getDataInfo().getLength()});
        data = new DataDoubleArray(dataInfo.getLength());
        dataInfo.addTag(tag);
        dataArray = data.getData();
    }
    
    /**
     * Returns the data source for volume scalings.
     */
    public IEtomicaDataSource getScalingDataSource() {
        return vDataSource;
    }
    
    public IData getData() {
        if (integrator == null) throw new IllegalStateException("must call setIntegrator before using meter");
        IBox box = integrator.getBox();
        for (int i=0; i<sim.getSpeciesCount(); i++) {
            myBox.setNMolecules(sim.getSpecies(i), box.getNMolecules(sim.getSpecies(i)));
        }
        IAtomList atoms = box.getLeafList();
        IAtomList myAtoms = myBox.getLeafList();
        for (int j=0; j<atoms.getAtomCount(); j++) {
            myAtoms.getAtom(j).getPosition().E(atoms.getAtom(j).getPosition());
        }
        inflater.setBox(myBox);
        energy.zeroSum();
        integrator.getPotentialMaster().calculate(box, iteratorDirective, energy);
        double uOld = energy.getSum();
        final double[] x = ((DataDoubleArray)vDataSource.getData()).getData();
        double mult = 1.0/space.D();
        
        IVector[] xNU = stuff.stuff();
        double vol0 = box.getBoundary().volume();
        IVectorMutable boxSize0 = space.makeVector();
        boxSize0.E(box.getBoundary().getBoxSize());
        ((PotentialMasterList)integrator.getPotentialMaster()).getNeighborManager(myBox).reset();
        IVectorMutable dr = space.makeVector();
        dr.Ev1Mv2(myAtoms.getAtom(0).getPosition(), myAtoms.getAtom(1).getPosition());
        box.getBoundary().nearestImage(dr);
        double r0 = Math.sqrt(dr.squared());
        
        for(int i=0; i<dataArray.length; i++) {
//System.out.println("scale "+Math.pow(x[i],mult));
            inflater.setScale(Math.pow(x[i],mult));
            inflater.actionPerformed();
            dr.Ev1Mv2(myAtoms.getAtom(0).getPosition(), myAtoms.getAtom(1).getPosition());
            box.getBoundary().nearestImage(dr);
            double r1 = Math.sqrt(dr.squared());
            
            for (int j=0; j<atoms.getAtomCount(); j++) {
                myAtoms.getAtom(j).getPosition().PEa1Tv1(vol0*(x[i]-1), xNU[j]);
            }
            if (i==0) {
//                System.out.println("MPBV "+atoms.getAtom(0).getPosition()+" => "+myAtoms.getAtom(0).getPosition());
//                System.out.println("  "+(x[i]-1)*vol0+" "+xNU[0]);
                dr.Ev1Mv2(myAtoms.getAtom(0).getPosition(), myAtoms.getAtom(1).getPosition());
                box.getBoundary().nearestImage(dr);
                double r2 = Math.sqrt(dr.squared());
//                System.out.print(String.format("%8.5f %8.5f %8.5f\n", r0, r1, r2));
//                System.out.println("    "+xNU[0]);
            }
            
            energy.zeroSum();
            integrator.getPotentialMaster().calculate(myBox, iteratorDirective, energy);
            double uNew = energy.getSum();
//            System.out.print(String.format("stuff % 5.4f  % 6.3f\n", x[i], (uNew-uOld)));
            dataArray[i] = Math.exp(-(uNew-uOld)/integrator.getTemperature()
                              + box.getMoleculeList().getMoleculeCount()*(x[i]-1));

            myBox.getBoundary().setBoxSize(boxSize0);
            for (int j=0; j<atoms.getAtomCount(); j++) {
                myAtoms.getAtom(j).getPosition().E(atoms.getAtom(j).getPosition());
            }
        }
//        System.out.println(Math.log(dataArray[0])/((x[0]-1)*vol0)+" "+Math.log(dataArray[1])/((x[1]-1)*vol0));

        return data;
    }

    private static final long serialVersionUID = 1L;
    protected ISimulation sim;
    protected IBox myBox;
    private DataDoubleArray data;
    private IEtomicaDataInfo dataInfo;
    private final DataTag tag;
    private double[] dataArray;
    private final BoxInflate inflater;
    private final IteratorDirective iteratorDirective;
    private final PotentialCalculationEnergySum energy = new PotentialCalculationEnergySum();
    private final ISpace space;
    private DataSourceUniform xDataSource;
    protected DataSourceExp vDataSource;
    private IntegratorBox integrator;
    protected IStuff stuff;
    public int count, countOverlap;
    
    /**
     * Transforms the scaling from linear (-s to +s) to exponential (exp(-s) to exp(+s))
     */
    protected static class DataSourceExp implements IEtomicaDataSource {
        public DataSourceExp(DataSourceUniform wrappedDataSource) {
            this.wrappedDataSource = wrappedDataSource;
            tag = new DataTag();
            int n = wrappedDataSource.getDataInfo().getLength();
            if (n % 2 == 1) n--;
            data = new DataDoubleArray(n);
            dataInfo = new DataInfoDoubleArray("something", Null.DIMENSION, new int[]{n});
        }
        
        public DataTag getTag() {
            return tag;
        }
        
        public IEtomicaDataInfo getDataInfo() {
            return dataInfo;
        }
        
        public void reset() {
            double[] x = ((DataDoubleArray)wrappedDataSource.getData()).getData();
            double[] scaledX = data.getData();
            int xIndex = 0;
            double dx = (wrappedDataSource.getXMax()-wrappedDataSource.getXMin())/wrappedDataSource.getNValues();
            for (int i=0; i<scaledX.length; i++, xIndex++) {
                // skip x=0 in the center
                if (Math.abs(x[xIndex]) < 0.25*dx) xIndex++;
                scaledX[i] = Math.exp(x[xIndex]);
            }
        }
        
        public IData getData() {
            return data;
        }
        
        protected final DataSourceUniform wrappedDataSource;
        protected final DataTag tag;
        protected DataInfoDoubleArray dataInfo;
        protected DataDoubleArray data;
    }
}
