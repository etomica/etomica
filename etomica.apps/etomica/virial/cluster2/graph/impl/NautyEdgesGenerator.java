package etomica.virial.cluster2.graph.impl;

import java.io.BufferedReader;
import java.io.IOException;

import etomica.virial.cluster2.graph.*;
import etomica.virial.cluster2.nauty.ProcessWrapper;

public class NautyEdgesGenerator extends AbstractEdgesGenerator {

  private static final String FLAG_NAUTY = "Nauty";
  private EdgesDecoder decoder;
  private ProcessWrapper nauty;
  private BufferedReader nautyReader;

  // ***********************
  // * PUBLIC CONSTRUCTORS
  // ***********************

  public NautyEdgesGenerator(GraphSet family, ProcessWrapper nautyProcess,
      EdgesFilter filter) {

    super(family, false, filter);
    nauty = nautyProcess;
    decoder = EdgesFactory.createSimpleEdgesDecoder(getNodeCount());
    computeTags();
  }

  public NautyEdgesGenerator(GraphSet family, ProcessWrapper nautyProcess) {

    this(family, nautyProcess, null);
  }

  // ***********************
  // * PROTECTED METHODS
  // ***********************

  protected void run() {

    try {
      nauty.run();
      nautyReader = new BufferedReader(nauty.getProcessOutput());
    } catch (IOException e) {
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

    return NautyEdgesGenerator.FLAG_NAUTY;
  }

  @Override
  protected void pop() {

    String line;
    // the enumeration is starting at this point, so run nauty
    if (!isStarted()) {
      run();
    }
    try {
      // first line: size of the automorphism group associated with the graph
      line = nautyReader.readLine();
      if (line == null || line.isEmpty()) {
        setTop(null);
        return;
      }
//      if (isStarted() && getTop() == null) {
//        return;
//      }
      // number of isomorphisms: N!/automorphism_group_size
      double coefficient = 1;
      int automorphismGroupSize = Integer.valueOf(line);
      if (automorphismGroupSize > 0) {
        for (int i = 1; i <= getNodeCount(); i++) {
          coefficient *= i;
        }
        coefficient /= automorphismGroupSize;
      }
      // second line: upper triangle encoding of the graph as a bit string; this
      // graph is the representative of its automorphism group
      setTop(EdgesFactory.createNautyEdges(nautyReader.readLine(), decoder,
          coefficient));
    } catch (IOException e) {
      setTop(null);
      return;
    }
  }
}