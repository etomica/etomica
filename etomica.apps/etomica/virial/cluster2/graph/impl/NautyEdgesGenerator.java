/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph.impl;

import java.io.BufferedReader;
import java.io.IOException;

import etomica.virial.cluster2.bitmap.BitmapFactory;
import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.EdgesFilter;
import etomica.virial.cluster2.graph.EdgesRepresentationFactory;
import etomica.virial.cluster2.graph.GraphFactory;
import etomica.virial.cluster2.nauty.ProcessWrapper;

public class NautyEdgesGenerator extends AbstractEdgesGenerator {

  private static final String TAG_NAUTY = "Nauty";
  private ProcessWrapper nauty;
  private BufferedReader nautyReader;
  private EdgesRepresentationFactory factory;

  public NautyEdgesGenerator(EdgesRepresentationFactory edgesFactory,
      ProcessWrapper nautyProcess, EdgesFilter filter) {

    super(false, filter);
    factory = edgesFactory;
    nauty = nautyProcess;
    computeTags();
  }

  public NautyEdgesGenerator(EdgesRepresentationFactory edgesFactory,
      ProcessWrapper nautyProcess) {

    this(edgesFactory, nautyProcess, null);
  }

  protected void run() {

    try {
      nauty.run();
      nautyReader = new BufferedReader(nauty.getProcessOutput());
    }
    catch (IOException e) {
      e.printStackTrace();
      throw new RuntimeException("Nauty IOException");
    }
  }

  @Override
  protected void computeTags() {

    getInternalTags().add(getTag());
    getInternalTags().addAll(nauty.getProcessInfo().getTags());
    if (getEdgesFilter() != null) {
      getInternalTags().addAll(getEdgesFilter().getTags());
    }
  }

  @Override
  protected String getTag() {

    return NautyEdgesGenerator.TAG_NAUTY;
  }

  @Override
  protected Edges push() {

    String line;
    // the enumeration is starting now, so run the modified nauty
    if (!isStarted()) {
      run();
    }
    try {
      // first line: size of the automorphism group associated with the graph
      line = nautyReader.readLine();
      if (line == null || line.trim().length() == 0) {
        return null;
      }
      // number of isomorphisms: N!/automorphism_group_size
      int coefficient = 1;
      int autoGroupSize = Integer.valueOf(line);
      if (autoGroupSize > 0) {
        for (int i = 1; i <= factory.getNodeCount(); i++) {
          coefficient *= i;
        }
        coefficient /= autoGroupSize;
      }
      // second line: encoding of the graph as a bit string; this
      // graph is the representative of its automorphism group
      return GraphFactory.nautyEdges(factory.getRepresentation(BitmapFactory
          .getBitmap(nautyReader.readLine())), GraphFactory.defaultCoefficient(coefficient));
    }
    catch (IOException e) {
      return null;
    }
  }
}