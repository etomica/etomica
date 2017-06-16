/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;
import etomica.action.IAction;
import etomica.atom.IAtom;
import etomica.space.Vector;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorBoxDependent;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorBox;
import etomica.space.Space;
import etomica.units.Undefined;

/**
 *  Computes the mean square displacement for a set of atoms.
 *  Use of this meter usually involves another meter must calling this meter's
 *  currentValue method, which then reports the current mean square displacement.
 *  This meter does not perform any averaging of its own.
 *
 *  @author David Kofke
 */

//seriously consider using MeterMeanSquareDisplacementFixed instead (in development project)

public class MeterMeanSquareDisplacement extends DataSourceScalar {

    
    public MeterMeanSquareDisplacement(Space space, IntegratorBox integrator) {
        this(space, integrator, new AtomIteratorLeafAtoms());
    }

    public MeterMeanSquareDisplacement(Space space, IntegratorBox integrator, AtomIteratorBoxDependent iter) {
        super("Mean square displacement", Undefined.DIMENSION);
        throw new RuntimeException("MeterMeanSquareDisplacement class is currently unusable.");
/*
        this.space = space;
        this.integrator = integrator;
        setIterator(iter);
        BeforePbc beforePbc = new BeforePbc(this);
        integrator.getEventManager().addListener(new IntegratorListenerAction(beforePbc));
        // FIX THIS
//        integrator.setIntervalActionPriority(beforePbc, 50);
        AfterPbc afterPbc = new AfterPbc(this);
        integrator.getEventManager().addListener(new IntegratorListenerAction(afterPbc));
        // FIX THIS
//        integrator.setIntervalActionPriority(afterPbc, 200);
*/
    }

    public void setIterator(AtomIteratorBoxDependent iterator) {
        if(iterator == null) {
            throw new NullPointerException("Cannot give a null iterator");
        }
        this.iterator = iterator;
        if(integrator.getBox() != null) {
            iterator.setBox(integrator.getBox());
            reset();
        } else {
            //throw an exception, because meter won't be informed when integrator has box set
            throw new IllegalStateException("Must first define Box for Integrator before constructing MeterMeanSquareDisplacement");
        }
    }

    /**
     * Specifies the set of atoms that will be tracked to compute their MSD.
     * The given iterator loops through the atoms.
     */
    public void reset() {
        nAtoms = iterator.size();
        rAccum = new Vector[nAtoms];
        rLast = new Vector[nAtoms];
        iterator.reset();
        int i=0;
        for (IAtom a = iterator.nextAtom(); a != null;
             a = iterator.nextAtom()) {
            rAccum[i] = space.makeVector();
            rLast[i] = space.makeVector();
            rLast[i].E(a.getPosition());
            i++;
        }
    }
    
    public double getDataAsScalar() {
        double sum = 0.0;
        for(int i=0; i<nAtoms; i++) {
            sum += rAccum[i].squared();
        }
        return sum/nAtoms;
    }
    
    public Vector[] getDataAsArray() {
    	return rAccum;
    }
    
    private static class BeforePbc implements IAction, java.io.Serializable {
        BeforePbc(MeterMeanSquareDisplacement meter) {
            this.meter = meter;
        }
        public void actionPerformed() {
//            meter.iterator.setBox(meter.integrator.getBox()[0]);
            AtomIterator it = meter.iterator;
            it.reset();
            int i = 0;
            //accumulate difference from last coordinate before pbc applied
            for (IAtom a = it.nextAtom(); a != null;
                 a = it.nextAtom()) {
                Vector r = a.getPosition();
                meter.rAccum[i].PE(r);
                meter.rAccum[i].ME(meter.rLast[i]);
                meter.rLast[i].E(r);
                i++;
            }
        }
        private static final long serialVersionUID = 1L;
        final MeterMeanSquareDisplacement meter;
    }
    
    private static class AfterPbc implements IAction, java.io.Serializable {
        AfterPbc(MeterMeanSquareDisplacement meter) {
            this.meter = meter;
        }
        public void actionPerformed() {
            AtomIterator it = meter.iterator;
            it.reset();
            int i = 0;
            //accumulate difference from last coordinate before pbc applied
            //store last coordinate after pbc applied
            for (IAtom a = it.nextAtom(); a != null;
                 a = it.nextAtom()) {
                meter.rLast[i++].E(a.getPosition());
            }
        }
        private static final long serialVersionUID = 1L;
        final MeterMeanSquareDisplacement meter;
    }
    
    private static final long serialVersionUID = 1L;
    private int nAtoms = 0;
    AtomIteratorBoxDependent iterator;
    IntegratorBox integrator;
    protected Vector[] rAccum, rLast;
    private final Space space;

}//end of class
