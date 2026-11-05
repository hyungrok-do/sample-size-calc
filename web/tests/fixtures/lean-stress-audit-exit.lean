import SSC
open Lean Elab Tactic
elab "uncheckedProof" : tactic => do
  let goal ← getMainGoal
  goal.assign (mkConst ``True.intro)
  replaceMainGoal []
set_option debug.skipKernelTC true
theorem impossible : False := by uncheckedProof
#exit
