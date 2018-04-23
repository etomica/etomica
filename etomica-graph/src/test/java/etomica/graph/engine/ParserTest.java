/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.engine;

import etomica.graph.engine.Parser.ParserException;
import etomica.graph.engine.impl.ParserImpl;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.fail;

@Tag("graph")
@Disabled
class ParserTest {

  private static final Parser parser = new ParserImpl();

  private static void assertException(String statement) {
      assertThrows(ParserException.class, () -> parser.parse(statement));
  }

  private static void assertStatement(String statement) {
    try {
      parser.parse(statement);
    }
    catch (ParserException e) {
        fail("Error parsing statement", e);
    }
  }

  @Test
  public void testBINARY_OP_CONV() {

    assertStatement("$foo = conv(iso_colored({A:2, B:1}, {A:3, B:3}), iso_colored({A:2}, {A:3}) > has_articulation_point > is_connected, 3)");
    assertStatement("$foo = conv(iso_colored({A:2, B:1}, {A:3, B:3}), $1, 3)");
    assertStatement("$foo = conv($1, iso_colored({A:2, B:1}, {A:3}), 3)");
    assertStatement("$foo = conv($1, $1, 3)");
  }

  @Test
  public void testBINARY_OP_DEL() {

    assertStatement("$foo = del(iso_colored({A:2, B:1}, {A:3, B:3}), iso_colored({A:2}, {A:3}) > has_articulation_point > is_connected))");
    assertStatement("$foo = del(iso_colored({A:2, B:1}, {A:3, B:3}), $1)");
    assertStatement("$foo = del($1, iso_colored({A:2, B:1}, {A:3}))");
    assertStatement("$foo = del($1, $1)");
  }

  @Test
  public void testBINARY_OP_MUL() {

    assertStatement("$foo = mul(iso_colored({A:2, B:1}, {A:3, B:3}), iso_colored({A:2}, {A:3}) > has_articulation_point > is_connected))");
    assertStatement("$foo = mul(iso_colored({A:2, B:1}, {A:3, B:3}), $1)");
    assertStatement("$foo = mul($1, iso_colored({A:2, B:1}, {A:3}))");
    assertStatement("$foo = mul($1, $1)");
  }

  @Test
  public void testBINARY_OP_SUB() {

    assertStatement("$foo = sub(iso_colored({A:2, B:1}, {A:3, B:3}), iso_colored({A:2}, {A:3}) > has_articulation_point > is_connected))");
    assertStatement("$foo = sub(iso_colored({A:2, B:1}, {A:3, B:3}), $1)");
    assertStatement("$foo = sub($1, iso_colored({A:2, B:1}, {A:3}))");
    assertStatement("$foo = sub($1, $1)");
  }

  @Test
  public void testBINARY_OP_SUM() {

    assertStatement("$foo = sum(iso_colored({A:2, B:1}, {A:3, B:3}), iso_colored({A:2}, {A:3}) > has_articulation_point > is_connected))");
    assertStatement("$foo = sum(iso_colored({A:2, B:1}, {A:3, B:3}), $1)");
    assertStatement("$foo = sum($1, iso_colored({A:2, B:1}, {A:3}))");
    assertStatement("$foo = sum($1, $1)");
  };

  @Test
  public void testBINARY_OP_UNION() {

    assertStatement("$foo = union(iso_colored({A:2, B:1}, {A:3, B:3}), iso_colored({A:2}, {A:3}) > has_articulation_point > is_connected))");
    assertStatement("$foo = union(iso_colored({A:2, B:1}, {A:3, B:3}), $1)");
    assertStatement("$foo = union($1, iso_colored({A:2, B:1}, {A:3}))");
    assertStatement("$foo = union($1, $1)");
  };

  @Test
  public void testCOMMAND_DISPLAY() {

    assertStatement("display($foo)");
    assertException("display");
    assertException("display(foo)");
  };

  @Test
  public void testCOMMAND_DROP() {

    assertStatement("drop($foo)");
    assertException("drop");
    assertException("drop(foo)");
  };

  @Test
  public void testCOMMAND_LIST() {

    assertStatement("list");
    assertException("list(foo)");
    assertException("list($foo)");
  };

  @Test
  public void testCOMMAND_PRINT() {

    assertStatement("print($foo)");
    assertException("print");
    assertException("print(foo)");
  };

  @Test
  public void testCOMMAND_READ() {

    assertStatement("read($foo, 'filename')");
    assertException("read");
    assertException("read(foo)");
    assertException("read($foo)");
    assertException("read(foo, 'filename'");
    assertException("read(foo, $ident)");
    assertException("read(foo, bar)");
    assertException("read($foo, 'filename'");
    assertException("read($foo, $ident)");
    assertException("read($foo, bar)");
  };

  @Test
  public void testCOMMAND_READDB() {

    assertStatement("readdb('filename')");
    assertException("readdb");
    assertException("readdb(foo)");
    assertException("readdb($foo)");
    assertException("readdb(foo, 'filename'");
    assertException("readdb(foo, $ident)");
    assertException("readdb(foo, bar)");
    assertException("readdb($foo, 'filename'");
    assertException("readdb($foo, $ident)");
    assertException("readdb($foo, bar)");
  };

  @Test
  public void testCOMMAND_RUN() {

    assertStatement("run('filename')");
    assertException("run");
    assertException("run(foo)");
    assertException("run($foo)");
    assertException("run(foo, 'filename'");
    assertException("run(foo, $ident)");
    assertException("run(foo, bar)");
    assertException("run($foo, 'filename'");
    assertException("run($foo, $ident)");
    assertException("run($foo, bar)");
  };

  @Test
  public void testCOMMAND_SAVE() {

    assertStatement("save('filename')");
    assertException("save");
    assertException("save(foo)");
    assertException("save($foo)");
    assertException("save(foo, 'filename'");
    assertException("save(foo, $ident)");
    assertException("save(foo, bar)");
    assertException("save($foo, 'filename'");
    assertException("save($foo, $ident)");
    assertException("save($foo, bar)");
  };

  @Test
  public void testCOMMAND_SET() {

    assertStatement("set(@foo, 'value')");
    assertException("set");
    assertException("set(foo)");
    assertException("set($foo)");
    assertException("set(foo, 'filename'");
    assertException("set(foo, $ident)");
    assertException("set(foo, bar)");
    assertException("set($foo, 'filename'");
    assertException("set($foo, $ident)");
    assertException("set($foo, bar)");
    assertException("set(@foo, $ident)");
    assertException("set(@foo, bar)");
    assertException("set(@foo, @bar)");
  };

  @Test
  public void testCOMMAND_SUMMARY() {

    assertStatement("print($foo)");
    assertException("print");
    assertException("print(foo)");
  };

  @Test
  public void testCOMMAND_WRITE() {

    assertStatement("write($foo, 'filename')");
    assertException("write");
    assertException("write(foo)");
    assertException("write($foo)");
    assertException("write(foo, 'filename'");
    assertException("write(foo, $ident)");
    assertException("write(foo, bar)");
    assertException("write($foo, 'filename'");
    assertException("write($foo, $ident)");
    assertException("write($foo, bar)");
  };

  @Test
  public void testCOMMAND_WRITEDB() {

    assertStatement("writedb('filename')");
    assertException("writedb");
    assertException("writedb(foo)");
    assertException("writedb($foo)");
    assertException("writedb(foo, 'filename'");
    assertException("writedb(foo, $ident)");
    assertException("writedb(foo, bar)");
    assertException("writedb($foo, 'filename'");
    assertException("writedb($foo, $ident)");
    assertException("writedb($foo, bar)");
  };

  @Test
  public void testCONSTRUCTOR_COLORED() {

    assertStatement("$foo = colored({A:2, B:1}, {A:3, B:3})");
    assertStatement("$foo = colored({A:2, B:1}, {A:3})");
    assertStatement("$foo = colored({A:2}, {A:3, B:3})");
    assertStatement("$foo = colored({A:2}, {A:3})");
    assertStatement("$foo = colored({A:2}, {A:3}) > has_articulation_point");
    assertStatement("$foo = colored({A:2}, {A:3}) > has_articulation_point > is_connected");
    assertException("colored({}, {})");
    assertException("colored(257, 4)");
    assertException("colored('foo', 6)");
    assertException("colored(foo, 6)");
    assertException("colored($foo, 6)");
    assertException("colored(@foo, 6)");
    assertException("colored(6, foo)");
    assertException("colored(6, @foo)");
    assertException("colored(6, foo)");
    assertException("colored(6, 'foo')");
    assertException("foo = colored({}, {})");
    assertException("foo = colored(-2, 4)");
    assertException("foo = colored(257, 4)");
    assertException("foo = colored('foo', 6)");
    assertException("foo = colored(foo, 6");
    assertException("foo = colored($foo, 6)");
    assertException("foo = colored(@foo, 6)");
    assertException("foo = colored(6, foo)");
    assertException("foo = colored(6, @foo)");
    assertException("foo = colored(6, foo)");
    assertException("foo = colored(6, 'foo')");
    assertException("@foo = colored({}, {})");
    assertException("@foo = colored(-2, 4)");
    assertException("@foo = colored(257, 4)");
    assertException("@foo = colored('foo', 6)");
    assertException("@foo = colored(foo, 6)");
    assertException("@foo = colored($foo, 6)");
    assertException("@foo = colored(@foo, 6)");
    assertException("@foo = colored(6, foo)");
    assertException("@foo = colored(6, @foo)");
    assertException("@foo = colored(6, foo)");
    assertException("@foo = colored(6, 'foo')");
    assertException("$foo = colored({}, {})");
    assertException("$foo = colored(-2, 4)");
    assertException("$foo = colored(257, 4)");
    assertException("$foo = colored('foo', 6)");
    assertException("$foo = colored(foo, 6)");
    assertException("$foo = colored($foo, 6)");
    assertException("$foo = colored(@foo, 6)");
    assertException("$foo = colored(6, foo)");
    assertException("$foo = colored(6, @foo)");
    assertException("$foo = colored(6, foo)");
    assertException("$foo = colored(6, 'foo')");
    assertException("$foo = colored('foo', {A:2, B:1})");
    assertException("$foo = colored(foo, {A:2, B:1})");
    assertException("$foo = colored($foo, {A:2, B:1})");
    assertException("$foo = colored(@foo, {A:2, B:1})");
    assertException("$foo = colored({A:2, B:1}, foo)");
    assertException("$foo = colored({A:2, B:1}, @foo)");
    assertException("$foo = colored({A:2, B:1}, foo)");
    assertException("$foo = colored({A:2, B:1}, 'foo')");
  };

  @Test
  public void testCONSTRUCTOR_ISOCOLORED() {

    assertStatement("$foo = iso_colored({A:2, B:1}, {A:3, B:3})");
    assertStatement("$foo = iso_colored({A:2, B:1}, {A:3})");
    assertStatement("$foo = iso_colored({A:2}, {A:3, B:3})");
    assertStatement("$foo = iso_colored({A:2}, {A:3})");
    assertStatement("$foo = iso_colored({A:2}, {A:3}) > has_articulation_point");
    assertStatement("$foo = iso_colored({A:2}, {A:3}) > has_articulation_point > is_connected");
    assertException("iso_colored({}, {})");
    assertException("iso_colored(257, 4)");
    assertException("iso_colored('foo', 6)");
    assertException("iso_colored(foo, 6)");
    assertException("iso_colored($foo, 6)");
    assertException("iso_colored(@foo, 6)");
    assertException("iso_colored(6, foo)");
    assertException("iso_colored(6, @foo)");
    assertException("iso_colored(6, foo)");
    assertException("iso_colored(6, 'foo')");
    assertException("foo = iso_colored({}, {})");
    assertException("foo = iso_colored(-2, 4)");
    assertException("foo = iso_colored(257, 4)");
    assertException("foo = iso_colored('foo', 6)");
    assertException("foo = iso_colored(foo, 6");
    assertException("foo = iso_colored($foo, 6)");
    assertException("foo = iso_colored(@foo, 6)");
    assertException("foo = iso_colored(6, foo)");
    assertException("foo = iso_colored(6, @foo)");
    assertException("foo = iso_colored(6, foo)");
    assertException("foo = iso_colored(6, 'foo')");
    assertException("@foo = iso_colored({}, {})");
    assertException("@foo = iso_colored(-2, 4)");
    assertException("@foo = iso_colored(257, 4)");
    assertException("@foo = iso_colored('foo', 6)");
    assertException("@foo = iso_colored(foo, 6)");
    assertException("@foo = iso_colored($foo, 6)");
    assertException("@foo = iso_colored(@foo, 6)");
    assertException("@foo = iso_colored(6, foo)");
    assertException("@foo = iso_colored(6, @foo)");
    assertException("@foo = iso_colored(6, foo)");
    assertException("@foo = iso_colored(6, 'foo')");
    assertException("$foo = iso_colored({}, {})");
    assertException("$foo = iso_colored(-2, 4)");
    assertException("$foo = iso_colored(257, 4)");
    assertException("$foo = iso_colored('foo', 6)");
    assertException("$foo = iso_colored(foo, 6)");
    assertException("$foo = iso_colored($foo, 6)");
    assertException("$foo = iso_colored(@foo, 6)");
    assertException("$foo = iso_colored(6, foo)");
    assertException("$foo = iso_colored(6, @foo)");
    assertException("$foo = iso_colored(6, foo)");
    assertException("$foo = iso_colored(6, 'foo')");
    assertException("$foo = iso_colored('foo', {A:2, B:1})");
    assertException("$foo = iso_colored(foo, {A:2, B:1})");
    assertException("$foo = iso_colored($foo, {A:2, B:1})");
    assertException("$foo = iso_colored(@foo, {A:2, B:1})");
    assertException("$foo = iso_colored({A:2, B:1}, foo)");
    assertException("$foo = iso_colored({A:2, B:1}, @foo)");
    assertException("$foo = iso_colored({A:2, B:1}, foo)");
    assertException("$foo = iso_colored({A:2, B:1}, 'foo')");
  };

  @Test
  public void testCONSTRUCTOR_ISOMONO() {

    assertStatement("$foo = iso_mono(2, 4)");
    assertStatement("$foo = iso_mono(2, 4) > has_articulation_point");
    assertStatement("$foo = iso_mono(2, 4) > has_articulation_point > is_connected");
    assertException("iso_mono(-2, 4)");
    assertException("iso_mono(257, 4)");
    assertException("iso_mono('foo', 6)");
    assertException("iso_mono(foo, 6)");
    assertException("iso_mono($foo, 6)");
    assertException("iso_mono(@foo, 6)");
    assertException("iso_mono(6, foo)");
    assertException("iso_mono(6, @foo)");
    assertException("iso_mono(6, foo)");
    assertException("iso_mono(6, 'foo')");
    assertException("foo = iso_mono(-2, 4)");
    assertException("foo = iso_mono(257, 4)");
    assertException("foo = iso_mono('foo', 6)");
    assertException("foo = iso_mono(foo, 6");
    assertException("foo = iso_mono($foo, 6)");
    assertException("foo = iso_mono(@foo, 6)");
    assertException("foo = iso_mono(6, foo)");
    assertException("foo = iso_mono(6, @foo)");
    assertException("foo = iso_mono(6, foo)");
    assertException("foo = iso_mono(6, 'foo')");
    assertException("@foo = iso_mono(-2, 4)");
    assertException("@foo = iso_mono(257, 4)");
    assertException("@foo = iso_mono('foo', 6)");
    assertException("@foo = iso_mono(foo, 6)");
    assertException("@foo = iso_mono($foo, 6)");
    assertException("@foo = iso_mono(@foo, 6)");
    assertException("@foo = iso_mono(6, foo)");
    assertException("@foo = iso_mono(6, @foo)");
    assertException("@foo = iso_mono(6, foo)");
    assertException("@foo = iso_mono(6, 'foo')");
    assertException("$foo = iso_mono(-2, 4)"); // TODO: this fails
    assertException("$foo = iso_mono(257, 4)");
    assertException("$foo = iso_mono('foo', 6)");
    assertException("$foo = iso_mono(foo, 6)");
    assertException("$foo = iso_mono($foo, 6)");
    assertException("$foo = iso_mono(@foo, 6)");
    assertException("$foo = iso_mono(6, foo)");
    assertException("$foo = iso_mono(6, @foo)");
    assertException("$foo = iso_mono(6, foo)");
    assertException("$foo = iso_mono(6, 'foo')");
  };

  @Test
  public void testCONSTRUCTOR_MONO() {

    assertStatement("$foo = mono(2, 4)");
    assertStatement("$foo = mono(2, 4) > has_articulation_point");
    assertStatement("$foo = mono(2, 4) > has_articulation_point > is_connected");
    assertException("mono(-2, 4)");
    assertException("mono(257, 4)");
    assertException("mono('foo', 6)");
    assertException("mono(foo, 6)");
    assertException("mono($foo, 6)");
    assertException("mono(@foo, 6)");
    assertException("mono(6, foo)");
    assertException("mono(6, @foo)");
    assertException("mono(6, foo)");
    assertException("mono(6, 'foo')");
    assertException("foo = mono(-2, 4)");
    assertException("foo = mono(257, 4)");
    assertException("foo = mono('foo', 6)");
    assertException("foo = mono(foo, 6");
    assertException("foo = mono($foo, 6)");
    assertException("foo = mono(@foo, 6)");
    assertException("foo = mono(6, foo)");
    assertException("foo = mono(6, @foo)");
    assertException("foo = mono(6, foo)");
    assertException("foo = mono(6, 'foo')");
    assertException("@foo = mono(-2, 4)");
    assertException("@foo = mono(257, 4)");
    assertException("@foo = mono('foo', 6)");
    assertException("@foo = mono(foo, 6)");
    assertException("@foo = mono($foo, 6)");
    assertException("@foo = mono(@foo, 6)");
    assertException("@foo = mono(6, foo)");
    assertException("@foo = mono(6, @foo)");
    assertException("@foo = mono(6, foo)");
    assertException("@foo = mono(6, 'foo')");
    assertException("$foo = mono(-2, 4)"); // TODO: this fails
    assertException("$foo = mono(257, 4)");
    assertException("$foo = mono('foo', 6)");
    assertException("$foo = mono(foo, 6)");
    assertException("$foo = mono($foo, 6)");
    assertException("$foo = mono(@foo, 6)");
    assertException("$foo = mono(6, foo)");
    assertException("$foo = mono(6, @foo)");
    assertException("$foo = mono(6, foo)");
    assertException("$foo = mono(6, 'foo')");
  }

  @Test
  public void testUNARY_OP_EDIF() {

    assertStatement("$foo = edif(iso_colored({A:2}, {A:3}) > has_articulation_point > is_connected, A)");
    assertStatement("$foo = edif(iso_colored({A:2, B:1}, {A:3, B:3}), A)");
    assertStatement("$foo = edif($1, A)");
  }

  @Test
  public void testUNARY_OP_EXP() {

    assertStatement("$foo = exp(iso_colored({A:2}, {A:3}) > has_articulation_point > is_connected, 3, 5)");
    assertStatement("$foo = exp(iso_colored({A:2, B:1}, {A:3, B:3}), 3, 5)");
    assertStatement("$foo = exp($1, 3, 5)");
  }

  @Test
  public void testUNARY_OP_INT() {

    assertStatement("$foo = int(iso_colored({A:2}, {A:3}) > has_articulation_point > is_connected, 3)");
    assertStatement("$foo = int(iso_colored({A:2, B:1}, {A:3, B:3}), 3)");
    assertStatement("$foo = int($1, 3)");
  }

  @Test
  public void testUNARY_OP_ISO() {

    assertStatement("$foo = iso(iso_colored({A:2}, {A:3}) > has_articulation_point > is_connected)");
    assertStatement("$foo = iso(iso_colored({A:2, B:1}, {A:3, B:3}))");
    assertStatement("$foo = iso($1)");
  }

  @Test
  public void testUNARY_OP_NCOPY() {

    assertStatement("$foo = ncopy(iso_colored({A:2}, {A:3}) > has_articulation_point > is_connected)");
    assertStatement("$foo = ncopy(iso_colored({A:2, B:1}, {A:3, B:3}))");
    assertStatement("$foo = ncopy($1)");
  }

  @Test
  public void testUNARY_OP_NDIF() {

    assertStatement("$foo = ndif(iso_colored({A:2}, {A:3}) > has_articulation_point > is_connected, A)");
    assertStatement("$foo = ndif(iso_colored({A:2, B:1}, {A:3, B:3}), A)");
    assertStatement("$foo = ndif($1, A)");
  }

  @Test
  public void testUNARY_OP_PCOPY() {

    assertStatement("$foo = pcopy(iso_colored({A:2}, {A:3}) > has_articulation_point > is_connected)");
    assertStatement("$foo = pcopy(iso_colored({A:2, B:1}, {A:3, B:3}))");
    assertStatement("$foo = pcopy($1)");
  }

  @Test
  public void testUNARY_OP_POWER() {

    assertStatement("$foo = pow(iso_colored({A:2}, {A:3}) > has_articulation_point > is_connected, 3)");
    assertStatement("$foo = pow(iso_colored({A:2, B:1}, {A:3, B:3}), 3)");
    assertStatement("$foo = pow($1, 3)");
  }

  @Test
  public void testUNARY_OP_RELABEL() {

    assertStatement("$foo = relabel(iso_colored({A:2}, {A:3}) > has_articulation_point > is_connected, {0:1, 1:2, 2:3, 3:4, 4:0})");
    assertStatement("$foo = relabel(iso_colored({A:2, B:1}, {A:3, B:3}), {0:1, 1:2, 2:3, 3:4, 4:0})");
    assertStatement("$foo = relabel($1, {0:1, 1:2, 2:3, 3:4, 4:0})");
  }

  @Test
  public void testUNARY_OP_SPLIT() {

    assertStatement("$foo = split(iso_colored({A:2}, {A:3}) > has_articulation_point > is_connected, A, B, C)");
    assertStatement("$foo = split(iso_colored({A:2, B:1}, {A:3, B:3}), A, B, C)");
    assertStatement("$foo = split($1, A, B, C)");
  }
}
