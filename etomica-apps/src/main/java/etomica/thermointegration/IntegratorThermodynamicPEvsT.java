/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.thermointegration;
import etomica.action.activity.ActivityIntegrate;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataPump;
import etomica.data.DataSourceIndependent;
import etomica.data.DataSourceUniform;
import etomica.data.DataTag;
import etomica.data.DataSourceUniform.LimitType;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataGroup;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.IntegratorBox;
import etomica.units.dimensions.Temperature;


public class IntegratorThermodynamicPEvsT extends IntegratorThermodynamic implements DataSourceIndependent {

    public IntegratorThermodynamicPEvsT(IntegratorBox subIntegrator) {
        super(subIntegrator);
        setTemperatureDataSource(new DataSourceUniform("temperature", Temperature.DIMENSION, 11, 1, 2, LimitType.INCLUSIVE, LimitType.INCLUSIVE));
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(subIntegrator.getPotentialMaster(), subIntegrator.getBox());
        accumulatorPE = new AccumulatorAverageCollapsing();
        DataPump dataPump = new DataPump(meterPE, accumulatorPE);
        dataPumps.add(dataPump);
    }

    public void reset() throws ConfigurationOverlapException {
        super.reset();
        iTemperatureStep = 0;

        double temperature = temperatureDataSource.getData().getValue(iTemperatureStep);
        // we know it's an IntegratorBox because we took that in the
        // constructor
        ((IntegratorBox)subIntegrator).setTemperature(temperature);
    }

    protected void advanceParameters() {
        data.getData()[iTemperatureStep] = ((DataGroup)accumulatorPE.getData()).getData(accumulatorPE.AVERAGE.index).getValue(0);
        System.out.println(temperatureDataSource.getData().getValue(iTemperatureStep)+" "+data.getData()[iTemperatureStep]);
        if (iTemperatureStep == temperatureDataSource.getNValues()-1) {
            // haha, Integrator will keep running.  it can't stop itself.  It's out of control!
            return;
        }

        iTemperatureStep++;
        double temperature = temperatureDataSource.getData().getValue(iTemperatureStep);
        // we know it's an IntegratorBox because we took that in the
        // constructor
        ((IntegratorBox)subIntegrator).setTemperature(temperature);
    }

    /**
     * Sets the temperature data source to the one given.  Each temperature
     * in the temperature data source will be visited by the integrator.
     */
    public void setTemperatureDataSource(DataSourceUniform newTemperatureDataSource) {
        temperatureDataSource = newTemperatureDataSource;
        data = new DataFunction(new int[]{newTemperatureDataSource.getNValues()});
    }

    /**
     * Returns the temperature data source.
     */
    public DataSourceUniform getTemperatureSource() {
        return temperatureDataSource;
    }

    public int getIndependentArrayDimension() {
        return 1;
    }

    public DataDoubleArray getIndependentData(int i) {
        return (DataDoubleArray)temperatureDataSource.getData();
    }

    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return (DataInfoDoubleArray)temperatureDataSource.getDataInfo();
    }

    public DataTag getIndependentTag() {
        return temperatureDataSource.getTag();
    }

    private static final long serialVersionUID = 1L;
    protected int iTemperatureStep;
    protected DataFunction data;
    protected DataSourceUniform temperatureDataSource;
    protected AccumulatorAverage accumulatorPE;
    
    public static void main(String[] args) {
        LjMc3D sim = new LjMc3D();
        
        System.out.println("equilibrating");
        ActivityIntegrate ai = (ActivityIntegrate)sim.getController().getAllActions()[0];
        ai.setMaxSteps(10000);
        sim.getController().actionPerformed();
        
        System.out.println("equilibration finished");

        sim.getController().removeAction(ai);
        IntegratorThermodynamicPEvsT integratorPET = new IntegratorThermodynamicPEvsT((IntegratorBox)sim.getIntegrator());
        integratorPET.setTemperatureDataSource(new DataSourceUniform("temperature", Temperature.DIMENSION, 
                                                                     11, 1, 2, LimitType.INCLUSIVE, LimitType.INCLUSIVE));
        ActivityIntegrate newAI = new ActivityIntegrate(integratorPET);
        newAI.setMaxSteps(220000);
        integratorPET.setNumSubsteps(20000);
        integratorPET.setNumEquilibrationSubsteps(10);
        sim.getController().addAction(newAI);
        sim.getController().actionPerformed();
    }
}
