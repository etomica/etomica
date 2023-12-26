/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.data;

import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.integrator.IntegratorMD;
import etomica.space.Vector;
import etomica.util.Statefull;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Writer;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * Listener that can store configurations.  For the purposes of MSD, it
 * will save all configurations that will be needed to compute MSD from.
 * With LOG2 storage, it will store configurations that are (more or less)
 * 1, 2, 4, 8, 16... steps ago.
 *
 * Configs are saved as part of the integrator's stepStarted event.
 */
public class ConfigurationStorage implements IntegratorListener, Statefull {

    public enum StorageType {LOG2, MSD, LINEAR}

    protected final Box box;
    protected Vector[][] configList;
    protected Vector[][] configVelList;
    protected long stepCount;
    protected final long[] savedSteps;
    protected final double[] savedTimes;
    protected StorageType storageType;
    protected final Vector[] dr;
    protected final Vector dri;
    protected final Set<ConfigurationStorageListener> listeners;
    protected boolean enabled;
    protected int interval, intervalCountdown;
    protected boolean doVel;
    protected double dt;

    public ConfigurationStorage(Box box, StorageType storageType) {
        this(box, storageType, 60, 1);
    }

    public ConfigurationStorage(Box box, StorageType storageType, int maxStored, int interval) {
        this.interval = interval;
        intervalCountdown = interval;
        this.box = box;
        this.storageType = storageType;
        savedSteps = new long[maxStored];
        savedTimes = new double[maxStored];
        for (int i = 0; i < savedSteps.length; i++) savedTimes[i] = savedSteps[i] = -1;
        dr = box.getSpace().makeVectorArray(box.getLeafList().size());
        dri = box.getSpace().makeVector();
        listeners = new HashSet<>();
        enabled = true;
        this.doVel = false;
    }

    public void setDoVelocity(boolean doVelocity) {
        doVel = doVelocity;
    }

    public void setSampleInterval(int newSampleInterval) {
        interval = intervalCountdown = newSampleInterval;
        reset();
    }

    public int getSampleInterval() {
        return interval;
    }

    public Box getBox() {
        return box;
    }

    public void setEnabled(boolean newEnabled) {
        enabled = newEnabled;
    }

    public boolean getEnabled() {
        return enabled;
    }

    public void reset() {
        for (int i = 0; i < savedSteps.length; i++) savedTimes[i] = savedSteps[i] = -1;
        stepCount = 0;
        intervalCountdown = interval;
    }

    public void addListener(ConfigurationStorageListener l) {
        listeners.add(l);
    }

    @Override
    public void integratorInitialized(IntegratorEvent e) {
    }

    @Override
    public void integratorStepStarted(IntegratorEvent e) {
        if (!enabled) return;
        intervalCountdown--;
        if (intervalCountdown > 0) return;
        intervalCountdown = interval;
        IAtomList atoms = box.getLeafList();
        int n = atoms.size();
        if (configList == null) {
            configList = new Vector[1][];
            configList[0] = box.getSpace().makeVectorArray(n);
            if (doVel) {
                configVelList = new Vector[1][];
                configVelList[0] = box.getSpace().makeVectorArray(n);
            }
        } else if (storageType == StorageType.LOG2) {
            // we need to copy existing configs forward to make way for our new config

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
                        if (doVel) {
                            configVelList = Arrays.copyOf(configVelList, j + 1);
                            configVelList[j] = box.getSpace().makeVectorArray(n);
                        }
                    }
                    for (int i = 0; i < n; i++) {
                        configList[j][i].E(configList[j - 1][i]);
                        if (doVel) {
                            configVelList[j][i].E(configVelList[j - 1][i]);
                        }
                    }
                }
            }
        } else if (storageType == StorageType.LINEAR) {
            if (configList.length < savedSteps.length) {
                Vector[][] newConfigList = new Vector[configList.length + 1][];
                Vector[][] newConfigVelList = doVel ? new Vector[configVelList.length + 1][] : null;
                for (int i = 0; i < configList.length; i++) {
                    newConfigList[i + 1] = configList[i];
                    if (doVel) newConfigVelList[i + 1] = configVelList[i];
                }
                configList = newConfigList;
                configList[0] = box.getSpace().makeVectorArray(n);
                if (doVel) {
                    configVelList = newConfigVelList;
                    configVelList[0] = box.getSpace().makeVectorArray(n);
                }
            } else {
                Vector[] tmp = configList[configList.length - 1];
                Vector[] tmpV = doVel ? configVelList[configVelList.length - 1] : null;
                for (int i = configList.length - 1; i > 0; i--) {
                    configList[i] = configList[i - 1];
                    if (doVel) configVelList[i] = configVelList[i - 1];
                }
                configList[0] = tmp;
                if (doVel) configVelList[0] = tmpV;
            }
            for (int i = configList.length - 1; i > 0; i--) { //no
                savedSteps[i] = savedSteps[i - 1];
                savedTimes[i] = savedTimes[i - 1];
            }
        }

        if (storageType == StorageType.MSD && stepCount == 1) {
            savedSteps[1] = savedSteps[0];
            savedTimes[1] = savedTimes[0];
            configList = Arrays.copyOf(configList, 2);
            configList[1] = configList[0];
            configList[0] = box.getSpace().makeVectorArray(n);
            if (doVel) {
                configVelList = Arrays.copyOf(configVelList, 2);
                configVelList[1] = configVelList[0];
                configVelList[0] = box.getSpace().makeVectorArray(n);
            }
        }

        // savedSteps isn't integrator's steps, it's just relative steps for our own purposes
        savedSteps[0] = stepCount;
        Integrator integrator = e.getIntegrator();
        savedTimes[0] = integrator instanceof IntegratorMD ? ((IntegratorMD) integrator).getCurrentTime() : integrator.getStepCount();
        if (dt==0 && stepCount>0) dt = savedTimes[0] - savedTimes[1];
        Vector boxDim = box.getBoundary().getBoxSize();
        for (int i = 0; i < n; i++) {
            Vector p = atoms.get(i).getPosition();
            Vector v = ((IAtomKinetic) atoms.get(i)).getVelocity();

            if (stepCount > 0) {//Sabry : PBC not needed for Vel

                // unwrap previously-stored configs to match current configs in terms
                // of computing displacements without worying about PBC
                dri.Ev1Mv2(p, configList[1][i]);
                dri.DE(boxDim);
                for (int k = 0; k < dri.getD(); k++) {
                    dri.setX(k, ((int) Math.round(dri.getX(k))) * boxDim.getX(k));
                }
                if (!dri.isZero()) {
                    for (int j = 1; j < configList.length; j++) {
                        configList[j][i].PE(dri);
                    }
                }
            }
            configList[0][i].E(p);
            if (doVel) configVelList[0][i].E(v);
        }
        for (ConfigurationStorageListener csl : listeners) {
            csl.newConfigruation();
        }
        if (storageType == StorageType.MSD) {
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
                    if (doVel) {
                        configVelList = Arrays.copyOf(configVelList, j + 2);
                        configVelList[j + 1] = configVelList[j];
                        configVelList[j] = box.getSpace().makeVectorArray(n);
                    }
                }
                savedSteps[j] = savedSteps[0];
                savedTimes[j] = savedTimes[0];
//                System.out.println("=> "+Arrays.toString(savedSteps));
                for (int i = 0; i < n; i++) {
                    configList[j][i].E(configList[0][i]);
                    if (doVel) configVelList[j][i].E(configVelList[0][i]);
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

    public Vector[] getSavedVel(int idx) {
        return configVelList[idx];
    }

    /**
     * @returns the time interval between consective configs.
     */
    public double getDeltaT() {
        return dt;
    }

    @Override
    public void integratorStepFinished(IntegratorEvent e) {
    }

    public interface ConfigurationStorageListener {
        void newConfigruation();
    }

    @Override
    public void saveState(Writer fw) throws IOException {
        fw.write(getClass().getName()+"\n");
        fw.write(""+stepCount+" "+intervalCountdown+" "+dt+"\n");
        for (int i=0; i<savedSteps.length && savedSteps[i] > -1; i++) {
            fw.write(""+savedSteps[i]+" "+savedTimes[i]+"\n");
            for (int j=0; j<configList[i].length; j++) {
                Vector p = configList[i][j];
                int D = p.getD();
                fw.write(""+p.getX(0));
                for (int k=1; k<D; k++) fw.write(" "+p.getX(k));
                if (doVel) {
                    Vector v = configVelList[i][j];
                    for (int k=0; k<D; k++) fw.write(" "+v.getX(k));
                }
                fw.write("\n");
            }
        }
        fw.write("\n");
    }

    @Override
    public void restoreState(BufferedReader br) throws IOException {
        if (!br.readLine().equals(getClass().getName())) {
            throw new RuntimeException("oops");
        }
        String[] bits = br.readLine().split(" ");
        stepCount = Long.parseLong(bits[0]);
        intervalCountdown = Integer.parseInt(bits[1]);
        dt = Double.parseDouble(bits[2]);
        int n = box.getLeafList().size();
        for (int i=0; i<savedSteps.length; i++) {
            String s = br.readLine();
            if (s.length() < 2) break;
            bits = s.split(" ");
            savedSteps[i] = Long.parseLong(bits[0]);
            savedTimes[i] = Double.parseDouble(bits[1]);
            if (configList == null) {
                configList = new Vector[0][];
                if (doVel) configVelList = new Vector[0][];
            }
            configList = Arrays.copyOf(configList, configList.length+1);
            configList[i] = box.getSpace().makeVectorArray(n);
            if (doVel) {
                configVelList = Arrays.copyOf(configVelList, configVelList.length+1);
                configVelList[i] = box.getSpace().makeVectorArray(n);
            }
            for (int j=0; j<configList[i].length; j++) {
                Vector p = configList[i][j];
                int D = p.getD();
                bits = br.readLine().split(" ");
                for (int k=0; k<D; k++) p.setX(k, Double.parseDouble(bits[k]));
                if (doVel) {
                    Vector v = configVelList[i][j];
                    for (int k=0; k<D; k++) v.setX(k, Double.parseDouble(bits[D+k]));
                }
            }
        }
    }
}