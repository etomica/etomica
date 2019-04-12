/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass2d;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.integrator.IntegratorMD;
import etomica.space.Vector;

import java.util.Arrays;

/**
 * Listener that can store configurations.  For the purposes of MSD, it
 * will save all configurations that will be needed to compute MSD from.
 * Otherwise, it will store configurations that are (more or less)
 * 1, 2, 4, 8, 16... steps ago.
 */
public class ConfigurationStorage implements IntegratorListener {
    protected final Box box;
    protected Vector[][] configList;
    protected long stepCount;
    protected long[] savedSteps = new long[60]; // can handle up to 2^60 steps
    protected double[] savedTimes = new double[60];
    protected boolean msd;
    protected final Vector[] dr;
    protected final Vector dri;

    public ConfigurationStorage(Box box, boolean storeForMSD) {
        this.box = box;
        this.msd = storeForMSD;
        for (int i = 0; i < savedSteps.length; i++) savedTimes[i] = savedSteps[i] = -1;
        dr = box.getSpace().makeVectorArray(box.getLeafList().size());
        dri = box.getSpace().makeVector();
    }

    public void reset() {
        for (int i = 0; i < savedSteps.length; i++) savedTimes[i] = savedSteps[i] = -1;
    }

    @Override
    public void integratorInitialized(IntegratorEvent e) {
    }

    @Override
    public void integratorStepStarted(IntegratorEvent e) {
        IAtomList atoms = box.getLeafList();
        int n = atoms.size();
        if (configList == null) {
            configList = new Vector[1][];
            configList[0] = box.getSpace().makeVectorArray(n);
        } else if (!msd) {
            // element j stores the config that is at least 1<<(j-1) steps ago (1 step ago, 2 steps ago, 4 steps ago)
            // at some point (step 8) all configList will store what is exacly 1<<(j-1) steps ago (0, -1, -2, -4, -8...)
            // the next step (9) will need to copy all elements forward (0, -1, -2, -3, -5, -9)
            // the next step (10) would be (0, -1, -2, -4, -6, -10)
            // then (11) (0, -1, -2, -3, -7, -11)
            // then (12) (0, -1, -2, -4, -8, -12)

            for (int j = configList.length; j > 0; j--) {
                if (1L << (j - 1) > stepCount) continue;
                // decide if we should copy j-1 into j
                if (stepCount % (1 << (j - 1)) == 0) {
                    // copy j-1 to j
                    savedSteps[j] = savedSteps[j - 1];
                    savedTimes[j] = savedTimes[j - 1];
                    if (configList.length <= j) {
                        configList = Arrays.copyOf(configList, j + 1);
                        configList[j] = box.getSpace().makeVectorArray(n);
                    }
                    for (int i = 0; i < n; i++) {
                        configList[j][i].E(configList[j - 1][i]);
                    }
                }
            }
        }

        if (msd && stepCount == 1) {
            savedSteps[1] = savedSteps[0];
            savedTimes[1] = savedTimes[0];
            configList = Arrays.copyOf(configList, 2);
            configList[1] = configList[0];
            configList[0] = box.getSpace().makeVectorArray(n);
        }

        // savedSteps isn't integrator's steps, it's just relative steps for our own purposes
        savedSteps[0] = stepCount;
        Integrator integrator = e.getIntegrator();
        savedTimes[0] = integrator instanceof IntegratorMD ? ((IntegratorMD) integrator).getCurrentTime() : integrator.getStepCount();
        Vector boxDim = box.getBoundary().getBoxSize();
        for (int i = 0; i < n; i++) {
            Vector p = atoms.get(i).getPosition();
            if (stepCount > 0) {
                dri.Ev1Mv2(p, configList[1][i]);
                dri.DE(boxDim);
                for (int k = 0; k < 2; k++) {
                    dri.setX(k, ((int) Math.round(dri.getX(k))) * boxDim.getX(k));
                }
                if (!dri.isZero()) {
                    for (int j = 1; j < configList.length; j++) {
                        configList[j][i].PE(dri);
                    }
                }
            }
            configList[0][i].E(p);
        }
        if (msd) {
            // copy our new config forward to whatever power of 2 we are
            // at step 32, we will save 5 copies of our current config.
            // our previous configs (like
//            System.out.println(stepCount+" "+Arrays.toString(savedSteps));
            for (int j = 1, d = 1; d < stepCount + 1 && stepCount % d == 0; j++, d *= 2) {
                if (configList.length <= j + 1) {
                    configList = Arrays.copyOf(configList, j + 2);
                    configList[j + 1] = configList[j];
                    savedSteps[j + 1] = savedSteps[j];
                    savedTimes[j + 1] = savedTimes[j];
                    configList[j] = box.getSpace().makeVectorArray(n);
                }
                savedSteps[j] = savedSteps[0];
                savedTimes[j] = savedTimes[0];
//                System.out.println("=> "+Arrays.toString(savedSteps));
                for (int i = 0; i < n; i++) {
                    configList[j][i].E(configList[0][i]);
                }
            }
        }

        stepCount++;
    }

    public int getLastConfigIndex() {
        if (savedSteps[0] == -1) return -1;
        for (int i = 0; i < savedSteps.length - 1; i++) {
            if (savedSteps[i + 1] == -1) return i;
        }
        return savedSteps.length - 1;
    }

    public double[] getSavedTimes() {
        return savedTimes;
    }

    public long[] getSavedSteps() {
        return savedSteps;
    }

    public Vector[] getSavedConfig(int idx) {
        return configList[idx];
    }

    @Override
    public void integratorStepFinished(IntegratorEvent e) {
    }
}
