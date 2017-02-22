/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.engine;

import java.util.Set;

import etomica.graph.engine.Parser.Assignment;
import etomica.graph.engine.Parser.BinaryOp;
import etomica.graph.engine.Parser.BinaryOpByte;
import etomica.graph.engine.Parser.Command;
import etomica.graph.engine.Parser.CommandPropertyValue;
import etomica.graph.engine.Parser.CommandValue;
import etomica.graph.engine.Parser.CommandVariable;
import etomica.graph.engine.Parser.CommandVariableValue;
import etomica.graph.engine.Parser.Constructor;
import etomica.graph.engine.Parser.ConstructorColored;
import etomica.graph.engine.Parser.ConstructorMono;
import etomica.graph.engine.Parser.Expression;
import etomica.graph.engine.Parser.ParserException;
import etomica.graph.engine.Parser.Statement;
import etomica.graph.engine.Parser.StrictExpression;
import etomica.graph.engine.Parser.UnaryOp;
import etomica.graph.engine.Parser.UnaryOpByte;
import etomica.graph.engine.Parser.UnaryOpByteMap;
import etomica.graph.engine.Parser.UnaryOpColor;
import etomica.graph.engine.Parser.UnaryOpThreeColor;
import etomica.graph.engine.Parser.UnaryOpTwoByte;
import etomica.graph.engine.Parser.Variable;
import etomica.graph.iterators.DefaultIterator;
import etomica.graph.iterators.IsomorphismPrefilteredPartitionedIterator;
import etomica.graph.iterators.IteratorToSet;
import etomica.graph.iterators.PartitionedIterator;
import etomica.graph.iterators.StoredIterator;
import etomica.graph.iterators.filters.IsomorphismFilter;
import etomica.graph.iterators.filters.PropertyFilter;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphIterator;
import etomica.graph.operations.Conv;
import etomica.graph.operations.ConvParameters;
import etomica.graph.operations.Delete;
import etomica.graph.operations.DifByEdge;
import etomica.graph.operations.DifByNode;
import etomica.graph.operations.DifParameters;
import etomica.graph.operations.Exp;
import etomica.graph.operations.ExpParameters;
import etomica.graph.operations.Int;
import etomica.graph.operations.IsoFree;
import etomica.graph.operations.Mul;
import etomica.graph.operations.Mul.MulParameters;
import etomica.graph.operations.NCopy;
import etomica.graph.operations.PCopy;
import etomica.graph.operations.Pow;
import etomica.graph.operations.PowParameters;
import etomica.graph.operations.Relabel;
import etomica.graph.operations.RelabelParameters;
import etomica.graph.operations.Split;
import etomica.graph.operations.SplitParameters;
import etomica.graph.operations.Sub;
import etomica.graph.operations.Sum;
import etomica.graph.operations.Union;
import etomica.graph.property.HasArticulationPair;
import etomica.graph.property.HasArticulationPoint;
import etomica.graph.property.HasNoRootEdge;
import etomica.graph.property.IsBiconnected;
import etomica.graph.property.IsConnected;
import etomica.graph.viewer.ClusterViewer;

public class Interpreter implements ConsoleReader {

  private Console console;
  private Environment environment;
  private Parser parser;
  private ConsoleContainer container;
  private long stop;
  private long start;

  private void notImplemented() {

    console.write("not implemented");
    console.writeLn();
  }

  private void echo() {

    console.writeLn();
  }

  private void echoLn(String value) {

    console.write(value);
    echo();
  }

  private void echo(String value) {

    console.write(value);
  }

  private void echoCommand(Command command) {

    console.write(command.getCommand());
  }

  // interpreter specific exception type
  public class InterpreterException extends Exception {

    private static final long serialVersionUID = -4211862381484158892L;

    public InterpreterException(String message) {

      super(message);
    }

    public InterpreterException(String message, Throwable cause) {

      super(message, cause);
    }
  }

  public Interpreter(ConsoleContainer container, Environment environment, Parser parser) {

    this.container = container;
    this.console = container.getConsole();
    this.environment = environment;
    this.parser = parser;
    console.setReader(this);
  }

  // the console sends individual statements (lines) to the listener
  public void read(String statement) {

    try {
      echo();
      Statement result = parser.parse(statement);
      if (result instanceof Command) {
        dispatch((Command) result);
      }
      else if (result instanceof Assignment) {
        dispatch((Assignment) result);
      }
    }
    catch (ParserException e) {
      echo(e.getMessage());
    }
    if (!Parser.COMMAND_CLEAR.equals(statement)) {
      echo();
    }
  }

  private void dispatch(Command cmd) {

    String command = cmd.getCommand();
    if (Parser.COMMAND_DISPLAY.equals(command)) {
      display((CommandVariable) cmd);
    }
    else if (Parser.COMMAND_DROP.equals(command)) {
      drop((CommandVariable) cmd);
    }
    else if (Parser.COMMAND_CLEAR.equals(command)) {
      clear(cmd);
    }
    else if (Parser.COMMAND_LIST.equals(command)) {
      list(cmd);
    }
    else if (Parser.COMMAND_QUIT.equals(command)) {
      quit(cmd);
    }
    else if (Parser.COMMAND_PRINT.equals(command)) {
      print((CommandVariable) cmd);
    }
    else if (Parser.COMMAND_READ.equals(command)) {
      read((CommandVariableValue) cmd);
    }
    else if (Parser.COMMAND_READDB.equals(command)) {
      readdb((CommandValue) cmd);
    }
    else if (Parser.COMMAND_RUN.equals(command)) {
      run((CommandValue) cmd);
    }
    else if (Parser.COMMAND_SAVE.equals(command)) {
      save((CommandValue) cmd);
    }
    else if (Parser.COMMAND_SET.equals(command)) {
      set((CommandPropertyValue) cmd);
    }
    else if (Parser.COMMAND_SUMMARY.equals(command)) {
      summary((CommandVariable) cmd);
    }
    else if (Parser.COMMAND_WRITE.equals(command)) {
      write((CommandVariableValue) cmd);
    }
    else if (Parser.COMMAND_WRITEDB.equals(command)) {
      writedb((CommandValue) cmd);
    }
  }

  private void dispatch(Assignment assignment) {

    Variable variable = assignment.getVariable();
    StrictExpression expression = assignment.getExpression();

    startTime();
    Set<Graph> value = buildGraphSet(expression);
    stopTime();
    environment.setVariable(variable.getName(), value);
    echo(String.format("%s graphs generated in %s.", getValueFmt(value.size()), getTimeFmt()));
  }

  private String getTimeFmt() {

    String[] units = { "ns", "µs", "ms", "sec", "min" };
    int index = 0;
    long elapsed = (stop - start);
    // µs
    if (elapsed > 1000) {
      elapsed = elapsed / 1000;
      index++;
      // ms
      if (elapsed > 1000) {
        elapsed = elapsed / 1000;
        index++;
        // sec
        if (elapsed > 1000) {
          elapsed = elapsed / 1000;
          index++;
          // min
          if (elapsed > 120) {
            elapsed = elapsed / 60;
            index++;
          }
        }
      }
    }
    return String.format("%d%s", elapsed, units[index]);
  }

  private String getValueFmt(int size) {

    String[] units = { "", "K", "M", "G" };
    int index = 0;
    long length = size;
    // K
    if (length > 1000) {
      length = length / 1000;
      index++;
      // M
      if (length > 1000) {
        length = length / 1000;
        index++;
        // G
        if (length > 1000) {
          length = length / 1000;
          index++;
        }
      }
    }
    return String.format("%d%s", length, units[index]);
  }

  private void stopTime() {

    this.stop = System.nanoTime();
  }

  private void startTime() {

    this.start = System.nanoTime();
  }

  private Set<Graph> buildGraphSet(Expression expression) {

    Set<Graph> result = null;
    if (expression instanceof Constructor) {
      result = constructor((Constructor) expression);
    }
    else if (expression instanceof UnaryOp) {
      result = unaryOp((UnaryOp) expression);
    }
    else if (expression instanceof BinaryOp) {
      result = binaryOp((BinaryOp) expression);
    }
    if (expression instanceof Variable) {
      result = environment.getVariableValue(((Variable) expression).getName());
    }
    return result;
  }

  private Set<Graph> binaryOp(BinaryOp expr) {

    if (expr instanceof BinaryOpByte) {
      BinaryOpByte op = (BinaryOpByte) expr;
      Conv conv = new Conv();
      ConvParameters params = new ConvParameters(op.getByte1(), new MulParameters((byte)100));
      return conv.apply(buildGraphSet(op.getExpression1()), buildGraphSet(op.getExpression2()), params);
    }
    else {
      if (expr.getOperation().equals(Parser.BINARY_OP_DEL)) {
        Delete del = new Delete();
        return del.apply(buildGraphSet(expr.getExpression1()), buildGraphSet(expr.getExpression2()), null);
      }
      else if (expr.getOperation().equals(Parser.BINARY_OP_MUL)) {
        Mul mul = new Mul();
        return mul.apply(buildGraphSet(expr.getExpression1()), buildGraphSet(expr.getExpression2()), null);
      }
      else if (expr.getOperation().equals(Parser.BINARY_OP_SUB)) {
        Sub sub = new Sub();
        return sub.apply(buildGraphSet(expr.getExpression1()), buildGraphSet(expr.getExpression2()), null);
      }
      else if (expr.getOperation().equals(Parser.BINARY_OP_SUM)) {
        Sum sum = new Sum();
        return sum.apply(buildGraphSet(expr.getExpression1()), buildGraphSet(expr.getExpression2()), null);
      }
      else if (expr.getOperation().equals(Parser.BINARY_OP_UNION)) {
        Union union = new Union();
        return union.apply(buildGraphSet(expr.getExpression1()), buildGraphSet(expr.getExpression2()), null);
      }
    }
    assert (false);
    return null;
  }

  private Set<Graph> unaryOp(UnaryOp expr) {

    if (expr.getOperation().equals(Parser.UNARY_OP_EDIF)) {
      UnaryOpColor op = (UnaryOpColor) expr;
      DifByEdge edif = new DifByEdge();
      DifParameters params = new DifParameters(op.getColor1());
      return edif.apply(buildGraphSet(expr.getExpression()), params);
    }
    else if (expr.getOperation().equals(Parser.UNARY_OP_EXP)) {
      UnaryOpTwoByte op = (UnaryOpTwoByte) expr;
      Exp exp = new Exp();
      ExpParameters params = new ExpParameters(op.getByte1(), op.getByte2());
      return exp.apply(buildGraphSet(expr.getExpression()), params);
    }
    else if (expr.getOperation().equals(Parser.UNARY_OP_INT)) {
      UnaryOpColor op = (UnaryOpColor) expr;
      Int integration = new Int();
      DifParameters params = new DifParameters(op.getColor1());
      return integration.apply(buildGraphSet(expr.getExpression()), params);
    }
    else if (expr.getOperation().equals(Parser.UNARY_OP_ISO)) {
      IsoFree iso = new IsoFree();
      return iso.apply(buildGraphSet(expr.getExpression()), null);
    }
    else if (expr.getOperation().equals(Parser.UNARY_OP_NCOPY)) {
      NCopy ncopy = new NCopy();
      return ncopy.apply(buildGraphSet(expr.getExpression()), null);
    }
    else if (expr.getOperation().equals(Parser.UNARY_OP_NDIF)) {
      UnaryOpColor op = (UnaryOpColor) expr;
      DifByNode ndif = new DifByNode();
      DifParameters params = new DifParameters(op.getColor1());
      return ndif.apply(buildGraphSet(expr.getExpression()), params);
    }
    else if (expr.getOperation().equals(Parser.UNARY_OP_PCOPY)) {
      PCopy pcopy = new PCopy();
      return pcopy.apply(buildGraphSet(expr.getExpression()), null);
    }
    else if (expr.getOperation().equals(Parser.UNARY_OP_POWER)) {
      UnaryOpByte op = (UnaryOpByte) expr;
      Pow pow = new Pow();
      PowParameters params = new PowParameters(op.getByte1());
      return pow.apply(buildGraphSet(expr.getExpression()), params);
    }
    else if (expr.getOperation().equals(Parser.UNARY_OP_RELABEL)) {
      UnaryOpByteMap op = (UnaryOpByteMap) expr;
      Relabel relabel = new Relabel();
      byte[] permutation = new byte[op.getColorMap().size()];
      for (Byte key : op.getColorMap().keySet()) {
        permutation[key] = op.getColorMap().get(key);
      }
      RelabelParameters params = new RelabelParameters(permutation);
      return relabel.apply(buildGraphSet(expr.getExpression()), params);
    }
    else if (expr.getOperation().equals(Parser.UNARY_OP_SPLIT)) {
      UnaryOpThreeColor op = (UnaryOpThreeColor) expr;
      Split split = new Split();
      SplitParameters params = new SplitParameters(op.getColor1(), op.getColor2(), op.getColor3());
      return split.apply(buildGraphSet(expr.getExpression()), params);
    }
    assert (false);
    return null;
  }

  private GraphIterator appendFilters(Constructor expr, GraphIterator gi) {

    GraphIterator result = new PropertyFilter(gi, new HasNoRootEdge());
    if (expr.hasFilter(Parser.FILTER_IS_CONNECTED)) {
      result = new PropertyFilter(gi, new IsConnected());
    }
    if (expr.hasFilter(Parser.FILTER_IS_BICONNECTED)) {
      result = new PropertyFilter(gi, new IsBiconnected());
    }
    if (expr.hasFilter(Parser.FILTER_HAS_ARTICULATION_POINT)) {
      result = new PropertyFilter(gi, new HasArticulationPoint());
    }
    if (expr.hasFilter(Parser.FILTER_HAS_ARTICULATION_PAIR)) {
      result = new PropertyFilter(gi, new HasArticulationPair());
    }
    return result;
  }

  private Set<Graph> constructor(Constructor expr) {

    GraphIterator gi = null;
    if (expr instanceof ConstructorMono) {
      ConstructorMono c = (ConstructorMono) expr;
      byte nodeCount = (byte) (c.getFieldNodes() + c.getRootNodes());
      if (c.isIsoFree() && nodeCount > 1 && nodeCount < 10 && c.getRootNodes() == 0) {
        gi = new StoredIterator(nodeCount);
        gi = appendFilters(c, gi);
      }
      else {
        gi = new DefaultIterator((byte) (c.getFieldNodes() + c.getRootNodes()), c.getRootNodes());
        gi = appendFilters(c, gi);
        if (c.isIsoFree()) {
          gi = new IsomorphismFilter(gi);
        }
      }
    }
    else if (expr instanceof ConstructorColored) {
      ConstructorColored c = (ConstructorColored) expr;
      if (c.isIsoFree()) {
        gi = new IsomorphismPrefilteredPartitionedIterator(c.getRootColorMap(), c.getFieldColorMap());
        gi = new IsomorphismFilter(appendFilters(c, gi));
      }
      else {
        gi = appendFilters(c, new PartitionedIterator(c.getRootColorMap(), c.getFieldColorMap()));
      }
    }
    IteratorToSet i2set = new IteratorToSet();
    return i2set.getSet(gi);
  }

  // executes the display command
  public void display(CommandVariable command) {

    String varName = command.getVariable().getName();
    Set<Graph> gs = environment.getVariableValue(varName);
    if (gs == null) {
      echo("Variable is not set.");
    }
    else {
      ClusterViewer.createView(varName, gs);
      echo(String.format("Displaying $%s.", varName));
    }
  }

  // executes the drop command
  public void drop(CommandVariable command) {

    String varName = command.getVariable().getName();
    environment.setVariable(varName, null);
    echo(String.format("Variable $%s was dropped.", varName));
  }

  // executes the clear command
  public void clear(Command command) {

    console.clear();
  }

  // executes the list command
  public void list(Command command) {

    int varCount = 0;
    int propCount = 0;
    for (String name : environment.getPropertyNames()) {
      echoLn(String.format("@%s = '%s'", name, environment.getPropertyValue(name)));
      propCount++;
    }
    for (String name : environment.getVariableNames()) {
      Set<Graph> gs = environment.getVariableValue(name);
      echoLn(String.format("$%s : graph set size is %s.", name, getValueFmt(gs.size())));
      varCount++;
    }
    echo();
    echo(String.format("%d properties and %d variables listed.", propCount, varCount));
  }

  // executes the print command
  public void print(CommandVariable command) {

    String varName = command.getVariable().getName();
    Set<Graph> gs = environment.getVariableValue(varName);
    if (gs == null) {
      echo("Variable is not set.");
    }
    else {
      console.updateBegin();
      try {
        int count = 0;
        echoLn(String.format("Printing graph set $%s.", varName));
        echo();
        for (Graph g : gs) {
          echoLn(g.toString());
          count++;
        }
        echo();
        echo(String.format("A total of %d graphs were printed.", count));
      }
      finally {
        console.updateDone();
      }
    }
  }

  // executes the quit command
  public void quit(Command command) {

    container.quit();
  }

  // executes the read command
  public void read(CommandVariableValue command) {

    echoCommand(command);
    notImplemented();
  }

  // executes the readdb command
  public void readdb(CommandValue command) {

    echoCommand(command);
    notImplemented();
  }

  // executes the run command
  public void run(CommandValue command) {

    echoCommand(command);
    notImplemented();
  }

  // executes the save command
  public void save(CommandValue command) {

    echoCommand(command);
    notImplemented();
  }

  // executes the set command
  public void set(CommandPropertyValue command) {

    String property = command.getProperty().getName().toUpperCase();
    String oldValue = environment.getPropertyValue(property);
    environment.setProperty(property, command.getValue().getValue());
    if (oldValue != null) {
      echo(String.format("Property @%s (old = '%s') was overwritten.", property, oldValue));
    }
    else {
      echo(String.format("Property @%s was set.", property));
    }
  }

  // executes the summary command
  public void summary(CommandVariable command) {

    Set<Graph> gs = environment.getVariableValue(command.getVariable().getName());
    if (gs == null) {
      echo("Variable is not set.");
    }
    else {
      echo(String.format("Graph set size is %s.", getValueFmt(gs.size())));
    }
  }

  // executes the write command
  public void write(CommandVariableValue command) {

    echoCommand(command);
    notImplemented();
  }

  // executes the writedb command
  public void writedb(CommandValue command) {

    echoCommand(command);
    notImplemented();
  }
}