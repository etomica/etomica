package etomica.graph.viewer;

import java.awt.Dimension;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import javax.swing.JFrame;

import etomica.graph.engine.Viewer;
import etomica.graph.model.Graph;

public class ClusterViewer extends JFrame implements Viewer {

  private static final long serialVersionUID = 4330455556112932184L;
  private static Map<String, Viewer> stock = new HashMap<String, Viewer>();

  private String name;

  private ClusterViewer(String name) {

    this.name = name;
    setTitle("Cluster Viewer: $" + name);
    setPreferredSize(new Dimension(800, 600));
    setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);

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

    TiledGraphs tg = new TiledGraphs(graphs, this);
    tg.draw();
  }
}