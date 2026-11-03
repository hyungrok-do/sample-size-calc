import Lean

open Lean

namespace DesignCompiler.Equation

/-- The rule registry is extended by ordinary Lean theorems in user modules. -/
initialize ruleExtension : SimplePersistentEnvExtension Name (Array Name) ←
  registerSimplePersistentEnvExtension {
    addEntryFn := Array.push
    addImportedFn := fun entries => entries.foldl (· ++ ·) #[]
  }

initialize
  registerBuiltinAttribute {
    name := `equation_rule
    descr := "An equality or equivalence used by the Lean equation compiler."
    add := fun declaration stx kind => do
      Attribute.Builtin.ensureNoArgs stx
      unless kind == .global do throwError "equation_rule must be global"
      let info ← getConstInfo declaration
      unless info.isTheorem do throwError "equation_rule requires a proved theorem"
      modifyEnv fun env => ruleExtension.addEntry env declaration
  }

end DesignCompiler.Equation
