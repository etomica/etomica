/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.mappedDensity.crystal;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.normalmode.CoordinateDefinition;
import etomica.space.Vector;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;

/**
 * Meter for calculating three dimensional MSD
 */
public class MeterMSD extends DataSourceScalar {


    protected final Box box;
    protected Vector rivector;
    protected double sum;

    /**
     * Vector describing the orientation of the profile.
     * For example, (1,0) is along the x-axis.
     */
    /**
     * Meter that defines the property being profiled.
     */
    protected CoordinateDefinition latticesite;

    /**
     * Default constructor sets profile along the y-axis, with 100 histogram points.
     */
    public MeterMSD(Box box, CoordinateDefinition latticesite) {
        super("MSD", new CompoundDimension(new Dimension[]{Length.DIMENSION}, new double[]{2}));
        this.box = box;
        this.latticesite = latticesite;
        this.rivector = box.getSpace().makeVector();
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    /**
     * Returns the profile for the current configuration.
     */
    @Override
    public double getDataAsScalar() {

        IAtomList atoms = box.getLeafList();
        double sum = 0;
        for (IAtom atom : atoms) {
            rivector.Ev1Mv2(atom.getPosition(), latticesite.getLatticePosition(atom));
            box.getBoundary().nearestImage(rivector);
             double r2i = rivector.squared();
      //       try {fw.write(""+r2i+"\n");}catch (IOException ex){throw new RuntimeException(ex);}
            sum += r2i;
         }

         return sum / atoms.size();
    }

    /**
     * @return Returns the box.
     */
    public Box getBox() {
        return box;
    }

}
