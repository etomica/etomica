package etomica.virial.cluster2.ui;

import javax.swing.JTree;
import javax.swing.border.EmptyBorder;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreeModel;

public class TreeView {

  public static JTree build() {

    JTree tree = new JTree(createSampleTreeModel());
// tree.putClientProperty(Options.TREE_LINE_STYLE_KEY,
// Options.TREE_LINE_STYLE_NONE_VALUE);
    tree.expandRow(3);
    tree.expandRow(2);
    tree.expandRow(1);
    tree.setBorder(new EmptyBorder(5, 5, 5, 5));
    return tree;
  }

  private static TreeModel createSampleTreeModel() {

    DefaultMutableTreeNode root = new DefaultMutableTreeNode("Projects");
    DefaultMutableTreeNode parent;
    //
    parent = new DefaultMutableTreeNode("Project 1");
    root.add(parent);
    parent.add(new DefaultMutableTreeNode("GS #1 [N=4, All Permutations]"));
    parent.add(new DefaultMutableTreeNode("GS #2 [N=4, Ree-Hoover]"));
    parent
        .add(new DefaultMutableTreeNode("GS #3 [N=4, No Articulation Pairs]"));
    //
    parent = new DefaultMutableTreeNode("Project 2");
    root.add(parent);
    parent.add(new DefaultMutableTreeNode("GS #1 [N=6, Biconnected]"));
    parent.add(new DefaultMutableTreeNode("GS #2 [N=6, Ree-Hoover]"));
    return new DefaultTreeModel(root);
  }
}