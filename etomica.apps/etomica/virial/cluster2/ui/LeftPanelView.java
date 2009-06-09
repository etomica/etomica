package etomica.virial.cluster2.ui;

import java.awt.BorderLayout;
import java.awt.Dimension;

import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.SwingConstants;
import javax.swing.border.EmptyBorder;

import com.jgoodies.looks.Options;
import com.jgoodies.uif_lite.component.Factory;
import com.jgoodies.uif_lite.panel.SimpleInternalFrame;

public class LeftPanelView {

  public static JComponent build() {

    JPanel panel = new JPanel(new BorderLayout());
    panel.setOpaque(false);
    panel.add(buildVerticalSplit());
    return panel;
  }

  private static JComponent buildVerticalSplit() {

    JSplitPane pane = Factory.createStrippedSplitPane(
        JSplitPane.VERTICAL_SPLIT, buildTopPanel(), buildBottomPanel(), 0.3f);
    pane.setOpaque(false);
    return pane;
  }

  private static JComponent buildTopPanel() {

    SimpleInternalFrame sif = new SimpleInternalFrame("Navigator");
    JScrollPane pane = Factory.createStrippedScrollPane(TreeView.build());
    pane.setBorder(new EmptyBorder(2, 2, 2, 2));
    sif.add(pane);
    sif.setPreferredSize(new Dimension(200, 400));
    return sif;
  }

  private static JComponent buildBottomPanel() {

    SimpleInternalFrame sif = new SimpleInternalFrame("Designer");
    JTabbedPane tabbedPane = new JTabbedPane(SwingConstants.TOP);
    tabbedPane.putClientProperty(Options.EMBEDDED_TABS_KEY, Boolean.TRUE);
    tabbedPane.addTab("Generation", buildPanel());
    tabbedPane.addTab("Manipulation", buildPanel());
    JScrollPane pane = Factory.createStrippedScrollPane(tabbedPane);
    pane.setBorder(new EmptyBorder(2, 2, 2, 2));
    sif.add(pane);
    sif.setPreferredSize(new Dimension(200, 200));
    return sif;
  }

  private static JComponent buildPanel() {

    JPanel panel = new JPanel();
// panel.setBackground(Color.WHITE);
    return panel;
  }
}