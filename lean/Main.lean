import DesignCompiler.Examples

open DesignCompiler Lean

def main (args : List String) : IO UInt32 := do
  let stdout ← IO.getStdout
  let stderr ← IO.getStderr
  match args with
  | [name] =>
    match Examples.getExample name with
    | some study => stdout.putStrLn (compile study).compress; return 0
    | none => stderr.putStrLn "Unknown example. Choose normal, influence, or ate."; return 2
  | _ =>
    stderr.putStrLn "usage: designc normal|influence|ate"
    return 2
