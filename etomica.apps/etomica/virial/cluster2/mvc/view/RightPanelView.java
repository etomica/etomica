/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc.view;

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

import etomica.virial.cluster2.graph.GraphSetFactory;
import etomica.virial.cluster2.graph.Nodes;

public class RightPanelView {

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

    SimpleInternalFrame sif = new SimpleInternalFrame("Cluster View");
    JTabbedPane tabbedPane = new JTabbedPane(SwingConstants.TOP);
    tabbedPane.setAutoscrolls(false);
    tabbedPane.putClientProperty(Options.EMBEDDED_TABS_KEY, Boolean.TRUE);
    SVGDraw d = new SVGDraw(GraphSetFactory.completeGraphSet(new char[] {
        Nodes.NODE_COLOR_1, Nodes.NODE_COLOR_1, Nodes.NODE_COLOR_1,
        Nodes.NODE_COLOR_1 }, new char[] { Nodes.NODE_COLOR_1 }, true));
    tabbedPane.addTab("SVG Cluster", d.getPanel());
    JScrollPane pane = Factory.createStrippedScrollPane(tabbedPane);
    pane.setBorder(new EmptyBorder(2, 2, 2, 2));
    sif.add(pane);
    sif.setPreferredSize(new Dimension(400, 450));
    return sif;
  }

  private static JComponent buildBottomPanel() {

    SimpleInternalFrame sif = new SimpleInternalFrame("Data View");
    JScrollPane pane = Factory.createStrippedScrollPane(TableView.build());
    pane.setBorder(new EmptyBorder(2, 2, 2, 2));
    sif.add(pane);
    sif.setPreferredSize(new Dimension(400, 250));
    return sif;
  }
}