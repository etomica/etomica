/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.types.DataFunction;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorBoxFasterer;

/**
 * DataProcessor that takes as input the average density and average Boltzmann
 * factor profiles and calculates the chemical potential profile as the sum of
 * the ideal gas contribution (kT ln(<rho>)) and the excess chemical
 * potential, (-kT ln(<exp(-beta U)>)).  The incoming data is expected to be
 * the Boltzmann factor, while the density is taken from a DataSource.
 *
 * @author Andrew Schultz
 */
public class DataProcessorChemicalPotential extends DataProcessor {

    public DataProcessorChemicalPotential(IntegratorBox integrator) {
        this.integrator = integrator;
    }

    public DataProcessorChemicalPotential(IntegratorBoxFasterer integrator) {
        this.integratorFasterer = integrator;
    }

    public void setDensityProfileDump(IDataSource newDensityProfileSource) {
        densityProfileSource = newDensityProfileSource;
    }

    public IDataSource getDenstiyProfileDump() {
        return densityProfileSource;
    }

    protected IData processData(IData inputData) {
        double[] oldY = ((DataFunction)inputData).getData();
        double[] newY = data.getData();
        IData densityData = densityProfileSource.getData();
        if (densityData == null) return null;
        double temp = integrator != null ? integrator.getTemperature() : integratorFasterer.getTemperature();
        for (int i=0; i<oldY.length; i++) {
            double density = densityData.getValue(i);
            if (density*oldY[i] == 0) {
                newY[i] = Double.NaN;
            }
            else {
                newY[i] = temp * Math.log(density/oldY[i]);
            }
        }
        return data;
    }

    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        data = new DataFunction(new int[]{inputDataInfo.getLength()});
        dataInfo = inputDataInfo.getFactory().makeDataInfo();
        dataInfo.addTag(tag);
        return dataInfo;
    }

    protected IDataSource densityProfileSource;
    protected DataFunction data;
    protected IntegratorBox integrator;
    protected IntegratorBoxFasterer integratorFasterer;
}
