/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.action.controller.Controller;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.data.ConfigurationStorage;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.space.Space;
import etomica.space.Vector;
import g3dsys.control.G3DSys;
import g3dsys.images.Ball;
import g3dsys.images.Bond;

import java.awt.*;

public class DisplayBoxCanvas3DGlass extends DisplayBoxCanvasG3DSys implements DisplayBoxCanvasGlass {

    protected ConfigurationStorage configStorage;
    protected int configIndex;
    protected final Vector dr;
    protected boolean drawDisplacement, flipDisplacement;
    protected final Ball[] oldBalls;
    protected final Bond[] bonds;

    public DisplayBoxCanvas3DGlass(DisplayBox _box, Space _space, Controller controller, ConfigurationStorage configStorage) {
        super(_box, _space, controller);
        dr = _space.makeVector();
        this.configStorage = configStorage;
        configIndex = 100;
        oldBalls = new Ball[_box.getBox().getLeafList().size()];
        bonds = new Bond[oldBalls.length];
    }

    @Override
    public void setConfigIndex(int idx) {
        configIndex = idx;
        repaint();
    }

    @Override
    public int getConfigIndex() {
        return configIndex;
    }

    public void setConfigStorage(ConfigurationStorage configStorage){
        this.configStorage = configStorage;
    }

    public ConfigurationStorage getConfigStorage(){return configStorage;}


    @Override
    public void setDrawDisplacement(boolean doDrawDisplacement) {
        if (this.drawDisplacement == doDrawDisplacement) return;
        this.drawDisplacement = doDrawDisplacement;
        DiameterHashGlass diameterHash = (DiameterHashGlass) displayBox.getDiameterHash();
        diameterHash.setFlipped(drawDisplacement && flipDisplacement);
        if (doDrawDisplacement) {
            IAtomList atoms = displayBox.getBox().getLeafList();
            for (int i = 0; i < atoms.size(); i++) {
                if (oldBalls[i] == null) {
                    Ball newBall = new Ball(gsys, G3DSys.getColix((displayBox.getColorScheme().getAtomColor(atoms.get(i)))), 0, 0, 0, 0);
                    gsys.addFig(newBall);
                    oldBalls[i] = newBall;
                    bonds[i] = new Bond(getG3DSys(), aam.getAgent(atoms.get(i)), oldBalls[i]);
                    gsys.addFig(bonds[i]);
                } else {
                    // check that our old and new atom are the same
                    if (aam.getAgent(atoms.get(i)) != bonds[i].getBall1()) {
                        gsys.removeFig(bonds[i]);
                        bonds[i] = new Bond(getG3DSys(), aam.getAgent(atoms.get(i)), oldBalls[i]);
                        gsys.addFig(bonds[i]);
                    }
                }
                oldBalls[i].setDrawable(true);
                bonds[i].setDrawable(true);
            }
        } else {
            for (int i = 0; i < bonds.length; i++) {
                oldBalls[i].setDrawable(false);
                bonds[i].setDrawable(false);
            }
        }
    }

    public void setFlipDisplacement(boolean flipDisplacement) {
        if (this.flipDisplacement == flipDisplacement) return;
        this.flipDisplacement = flipDisplacement;
        DiameterHashGlass diameterHash = (DiameterHashGlass) displayBox.getDiameterHash();
        diameterHash.setFlipped(drawDisplacement && flipDisplacement);
        if (drawDisplacement) {
            IAtomList atoms = displayBox.getBox().getLeafList();
            for (int i = 0; i < atoms.size(); i++) {
                if (flipDisplacement) {
                    oldBalls[i].setD((float) diameterHash.getActualDiameter(atoms.get(i)));
                    aam.getAgent(atoms.get(i)).setD(0);
                } else {
                    aam.getAgent(atoms.get(i)).setD((float) diameterHash.getActualDiameter(atoms.get(i)));
                    oldBalls[i].setD(0);
                }
            }
        }
    }

    public boolean getFlipDisplacement() {
        return flipDisplacement;
    }

    @Override
    public boolean getDrawDisplacement() {
        return drawDisplacement;
    }

    public void doPaint(Graphics g) {
        if (drawDisplacement) {
            int idx = configIndex;
            int lastIndex = configStorage.getLastConfigIndex();
            if (idx > lastIndex) idx = lastIndex;
            if (idx >= 1) {
                double[] coords = new double[3];
                IAtomList atoms = displayBox.getBox().getLeafList();
                Vector[] oldPositions = configStorage.getSavedConfig(idx);
                AtomTestDeviation atomTest = (AtomTestDeviation) displayBox.getAtomTestDoDisplay();
                for (int i = 0; i < atoms.size(); i++) {
                    IAtom a = atoms.get(i);
                    Vector rOld = oldPositions[a.getLeafIndex()];
                    rOld.assignTo(coords);
                    oldBalls[i].setColor(G3DSys.getColix(displayBox.getColorScheme().getAtomColor(a)));
                    oldBalls[i].setX((float) coords[0]);
                    oldBalls[i].setY((float) coords[1]);
                    oldBalls[i].setZ((float) coords[2]);
                    oldBalls[i].setDrawable(atomTest == null || atomTest.test(a));
                }
            }
        }
        super.doPaint(g);
    }

}
