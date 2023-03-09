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
import etomica.util.voro.*;
import g3dsys.control.G3DSys;
import g3dsys.images.Ball;
import g3dsys.images.Bond;
import g3dsys.images.Line;
import org.jmol.util.Point3f;

import java.awt.*;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

public class DisplayBoxCanvas3DGlass extends DisplayBoxCanvasG3DSys implements DisplayBoxCanvasGlass {

    protected ConfigurationStorage configStorage;
    protected int configIndex;
    protected final Vector dr;
    protected boolean drawDisplacement, flipDisplacement, showVoronoiCells;
    protected final Ball[] oldBalls;
    protected final Bond[] bonds;
    protected final ArrayList<Line> voronoiEdges;
    protected double[] radii;

    public DisplayBoxCanvas3DGlass(DisplayBox _box, Space _space, Controller controller, ConfigurationStorage configStorage) {
        super(_box, _space, controller);
        dr = _space.makeVector();
        this.configStorage = configStorage;
        configIndex = 100;
        oldBalls = new Ball[_box.getBox().getLeafList().size()];
        bonds = new Bond[oldBalls.length];
        voronoiEdges = new ArrayList<>();
    }

    public void setVoronoiRadii(double[] r) {
        radii = r;
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

    public void setShowVoronoiCells(boolean b) {
        showVoronoiCells = b;
    }

    public boolean getFlipDisplacement() {
        return flipDisplacement;
    }

    @Override
    public boolean getDrawDisplacement() {
        return drawDisplacement;
    }

    public boolean getShowVoronoiCells() {return showVoronoiCells;}

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
        for (Line l : voronoiEdges) {
            gsys.removeFig(l);
        }
        voronoiEdges.clear();
        if (showVoronoiCells) {
            Vector b = displayBox.getBox().getBoundary().getBoxSize();
            double bx = 0.5*b.getX(0);
            double by = 0.5*b.getX(1);
            double bz = 0.5*b.getX(2);

            PreContainerPoly pconp=new PreContainerPoly(-bx,bx,-by,by,-bz,bz,true,true,true);
            for (IAtom a : displayBox.getBox().getLeafList()) {
                Vector r = a.getPosition();
                pconp.put(a.getLeafIndex(), r.getX(0), r.getX(1), r.getX(2), radii[a.getType().getIndex()]);
            }

            int[] nxout = new int[1];
            int[] nyout = new int[1];
            int[] nzout = new int[1];
            pconp.guess_optimal(nxout,nyout,nzout);
            int nx = nxout[0];
            int ny = nyout[0];
            int nz = nzout[0];
            short cWhite = G3DSys.getColix(Color.WHITE);

            ContainerPoly con = new ContainerPoly(-bx,bx,-by,by,-bz,bz,nx,ny,nz,true,true,true, 8);
            pconp.setup(con);
            CLoopAll vl = new CLoopAll(con);
            int[] pid = new int[1];
            Set<VoronoiCellBase.Edge> allEdges = new HashSet<VoronoiCellBase.Edge>();
            if(vl.start()) {
                VoronoiCell c = new VoronoiCell(con);
                do {
                    if(con.compute_cell(c,vl)) {
                        double[] x = new double[1];
                        double[] y = new double[1];
                        double[] z = new double[1];
                        double[] r = new double[1];
                        vl.pos(pid,x,y,z,r);
                        java.util.List<VoronoiCellBase.Edge> edges = c.edges(x[0], y[0], z[0]);
                        for (VoronoiCellBase.Edge e : edges) {
                            if (allEdges.contains(e)) continue;
                            allEdges.add(e);
                            Point3f p1 = Point3f.new3((float)e.x1, (float)e.y1, (float)e.z1);
                            Point3f p2 = Point3f.new3((float)e.x2, (float)e.y2, (float)e.z2);
                            Line l = new Line(gsys, cWhite, p1, p2);
                            gsys.addFig(l);
                            voronoiEdges.add(l);
                        }
                    }
                } while(vl.inc());
            }
        }
        super.doPaint(g);
    }

}
