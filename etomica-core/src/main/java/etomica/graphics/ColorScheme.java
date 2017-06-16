/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;
import java.awt.Color;

import etomica.action.IAction;
import etomica.action.SimulationRestart;
import etomica.atom.IAtom;
import etomica.modifier.ModifierGeneral;
import etomica.nbr.list.PotentialMasterList;

/**
 * Class that defines the algorithm used to determine atoms colors when drawn to DisplayBox.
 * The atomColor method is called just before the atom is drawn to set the graphics color.
 *
 * @author David Kofke
 */
 
public abstract class ColorScheme {

    protected Color defaultColor;
    
    public ColorScheme() {
        this(DEFAULT_ATOM_COLOR);
    }
    public ColorScheme(Color color) {
        defaultColor = color;
    }

    public abstract Color getAtomColor(IAtom a);
    
    public final void setDefaultColor(Color c) {defaultColor = c;}
    public final Color getDefaultColor() {return defaultColor;}

    public final static Color DEFAULT_ATOM_COLOR = Color.RED;
    
    /**
     * Colors all atoms with baseColor.
     */
    public static class Simple extends ColorScheme {
        public Simple() {super();}
        public Simple(java.awt.Color color) {super(color);}
        public Color getAtomColor(IAtom a) {return defaultColor;}
    }
    
    //for color scheme demo
    public static void main(String args[]) {
      final String APP_NAME = "Color Scheme";

      etomica.space.Space sp = etomica.space3d.Space3D.getInstance();
      final etomica.simulation.prototypes.HSMD3D sim = new etomica.simulation.prototypes.HSMD3D(sp);
      final SimulationGraphic simGraphic = new SimulationGraphic(sim, APP_NAME, sp, sim.getController());

      IAction repaintAction = simGraphic.getPaintAction(sim.box);

      DeviceNSelector nSelector = new DeviceNSelector(sim.getController());
      nSelector.setResetAction(new SimulationRestart(sim, sp, sim.getController()));
      nSelector.setSpecies(sim.species);
      nSelector.setBox(sim.box);
      nSelector.setPostAction(repaintAction);
      simGraphic.add(nSelector);


      simGraphic.getController().getReinitButton().setPostAction(repaintAction);

      final ColorSchemeByType ct = new ColorSchemeByType(sim);
      final ColorSchemeTemperature ctemp = new ColorSchemeTemperature(0,5);
      final ColorSchemeColliders ccld = new ColorSchemeColliders(sim.integrator);
      final ColorSchemeNeighbor nghb = new ColorSchemeNeighbor(sim, (PotentialMasterList)sim.potentialMaster, sim.box);
      nghb.setAtom(sim.box.getLeafList().getAtom(0));
      final ColorSchemeRandom rand = new ColorSchemeRandom(sim.box, sim.getRandom());
      final ColorSchemeCell cell = new ColorSchemeCell((PotentialMasterList)sim.potentialMaster,sim.getRandom(),sim.box);
      cell.setLattice(((PotentialMasterList)sim.potentialMaster).getNbrCellManager(sim.box).getLattice());
      
      IAction act = new IAction() {
        public void actionPerformed() {
          DisplayBox dp = (DisplayBox)simGraphic.displayList().getFirst();
          ct.setColor(sim.species.getLeafType(), 
              new java.awt.Color(sim.getRandom().nextInt(256),
                  sim.getRandom().nextInt(256),
                  sim.getRandom().nextInt(256))
              );
          dp.setColorScheme(ct);
        }
      };
      IAction act2 = new IAction() {
        public void actionPerformed() {
          DisplayBox dp = (DisplayBox)simGraphic.displayList().getFirst();
          dp.setColorScheme(ctemp);
        }
      };
      IAction act3 = new IAction() {
        public void actionPerformed() {
          DisplayBox dp = (DisplayBox)simGraphic.displayList().getFirst();
          dp.setColorScheme(ccld);
        }
      };
      IAction act4 = new IAction() {
        public void actionPerformed() {
          DisplayBox dp = (DisplayBox)simGraphic.displayList().getFirst();
          dp.setColorScheme(nghb);
        }
      };
      IAction act5 = new IAction() {
        public void actionPerformed() {
          DisplayBox dp = (DisplayBox)simGraphic.displayList().getFirst();
          dp.setColorScheme(rand);
        }
      };
      IAction act6 = new IAction() {
        public void actionPerformed() {
          ct.setColor(sim.species.getLeafType(), java.awt.Color.red);
          DisplayBox dp = (DisplayBox)simGraphic.displayList().getFirst();
          dp.setColorScheme(ct);
        }
      };
      IAction act7 = new IAction() {
        public void actionPerformed() {
          DisplayBox dp = (DisplayBox)simGraphic.displayList().getFirst();
          dp.setColorScheme(cell);
        }
      };

      DeviceButton colorer = new DeviceButton(sim.getController(),act);
      colorer.setLabel("Global Random");
      colorer.setPostAction(repaintAction);
      DeviceButton tempcolorer = new DeviceButton(sim.getController(),act2);
      tempcolorer.setLabel("By Temperature (0-5)");
      DeviceButton colliders = new DeviceButton(sim.getController(),act3);
      colliders.setLabel("Colliders");
      DeviceButton neighbors = new DeviceButton(sim.getController(),act4);
      neighbors.setLabel("Neighbors");
      neighbors.setPostAction(repaintAction);
      DeviceButton randomcol = new DeviceButton(sim.getController(),act5);
      randomcol.setLabel("Unique Random");
      randomcol.setPostAction(repaintAction);
      DeviceButton def = new DeviceButton(sim.getController(),act6);
      def.setLabel("Default Red");
      def.setPostAction(repaintAction);
      DeviceButton cellbtn = new DeviceButton(sim.getController(),act7);
      cellbtn.setLabel("Cell");
      
      DeviceSlider slabslide = new DeviceSlider(sim.getController());
      slabslide.setMinimum(0);
      slabslide.setMaximum(100);
      final DisplayBoxCanvasG3DSys dpg3d = (DisplayBoxCanvasG3DSys)simGraphic.getDisplayBox(sim.box).canvas;
      ModifierGeneral mg = new ModifierGeneral(dpg3d,"slab");
      slabslide.setModifier(mg);
      slabslide.setLabel("Slab");
      
      DeviceSlider depthslide = new DeviceSlider(sim.getController());
      depthslide.setMinimum(0);
      depthslide.setMaximum(100);
      ModifierGeneral mg2 = new ModifierGeneral(dpg3d,"depth");
      depthslide.setModifier(mg2);
      depthslide.setLabel("Depth");
      
      simGraphic.add(def);
      simGraphic.add(colorer);
      simGraphic.add(tempcolorer);
      simGraphic.add(colliders);
      simGraphic.add(neighbors);
      simGraphic.add(randomcol);
      simGraphic.add(cellbtn);
      simGraphic.add(slabslide);
      simGraphic.add(depthslide);
      
      
//      d.get
//      Action setslab = new Action() {
//
//        public void actionPerformed() {
//          //d.setS
//        }
//      };
      
      simGraphic.makeAndDisplayFrame(APP_NAME);
      
      //ColorSchemeTemperature colorScheme = new ColorSchemeTemperature(0.0f,5.0f);

      ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());
      colorScheme.setColor(sim.species.getLeafType(), java.awt.Color.red);

      DisplayBox dp = (DisplayBox)simGraphic.displayList().getFirst();
      dp.setColorScheme(colorScheme);
    }
}
