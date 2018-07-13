/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.viewer;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.UIManager;
import javax.swing.border.EmptyBorder;

import etomica.graph.engine.Viewer;
import etomica.graph.model.Graph;

public class ClusterViewer extends JFrame implements Viewer {

  private static final long serialVersionUID = 4330455556112932184L;
  private static Map<String, Viewer> stock = new HashMap<String, Viewer>();

  private String name;
  private JPanel mainPanel;

  static {
    // work around Java 1.7 problem.
    // getFocusedWindow seems to always be null and so
    // ToolTipManager.showTipWindow() displays no tooltips.
    UIManager.put("ToolTipManager.enableToolTipMode", "allWindows");
  }

  private ClusterViewer(String name) {

    this.name = name;
    setTitle("Cluster Viewer: $" + name);
    setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
    mainPanel = new JPanel();
    mainPanel.setPreferredSize(new Dimension(800, 600));
    mainPanel.setBorder(new EmptyBorder(10, 10, 10, 10));
    mainPanel.setLayout(new BorderLayout(0, 0));
    add(mainPanel);
    pack();
  }

  public static void createView(String name, Set<Graph> graphs) {

    Viewer v = stock.get(name);
    if (v == null) {
      v = new ClusterViewer(name);
      stock.put(name, v);
    }
    v.update(graphs);
    v.open();
  }

  public void open() {

    setVisible(true);
  }

  public void close() {

    setVisible(false);
  }

  public void setVisible(boolean b) {

    super.setVisible(b);
    if (!b) {
      stock.remove(name);
      dispose();
    }
  }

  public void update(Set<Graph> graphs) {

//    GraphMap gm = new GraphMap(graphs, mainPanel);
//    gm.draw();
  }
}
