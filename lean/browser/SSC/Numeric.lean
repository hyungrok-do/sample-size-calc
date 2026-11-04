import Init

/-!
Reusable floating-point numerical mathematics. These algorithms are executable
approximations, not kernel-certified error bounds. Failure to converge is explicit.
There are no study names, analysis selectors, or sample-size formulas here.
-/
namespace SSC.Numeric

abbrev Result := Except String

def nf (n : Nat) : Float := n.toFloat
def clamp01 (x : Float) : Float := max 0 (min 1 x)
def finite (x : Float) : Result Float :=
  if x.isFinite then .ok x else .error "Non-finite numerical result."

private def lanczos (x : Float) : Float := Id.run do
  let cs : Array Float := #[676.5203681218851, -1259.1392167224028,
    771.3234287776531, -176.6150291621406, 12.507343278686905,
    -0.13857109526572012, 9.984369578019572e-6, 1.5056327351493116e-7]
  let z := x - 1
  let mut a := 0.99999999999980993
  for i in [:cs.size] do a := a + cs[i]! / (z + nf (i+1))
  let t := z + 7.5
  return 0.91893853320467274178 + (z + 0.5) * t.log - t + a.log

def logGamma (x : Float) : Result Float := do
  if !x.isFinite || x <= 0 then throw "Gamma shape must be finite and positive."
  finite (if x < 0.5 then lanczos (x+1) - x.log else lanczos x)

/-- Regularized incomplete gamma, series / modified Lentz fraction. -/
def gammaCDF (a x : Float) : Result Float := do
  if !a.isFinite || a <= 0 || !x.isFinite then throw "Invalid incomplete-gamma arguments."
  if x <= 0 then return 0
  if a > 100000 then throw "Incomplete-gamma shape exceeds the numerical limit (100000)."
  let factor := (a * x.log - x - (← logGamma a)).exp
  if x < a + 1 then
    let mut term := 1 / a
    let mut total := term
    for i in [1:10001] do
      term := term * x / (a + nf i)
      total := total + term
      if term.abs <= total.abs * 2e-15 then return clamp01 (← finite (factor * total))
  else
    let mut b := x + 1 - a
    let mut c := 1e300
    let mut d := 1 / b
    let mut h := d
    for i in [1:10001] do
      let an := -(nf i) * (nf i - a)
      b := b + 2
      d := an * d + b
      if d.abs < 1e-300 then d := 1e-300
      c := b + an / c
      if c.abs < 1e-300 then c := 1e-300
      d := 1 / d
      let delta := d * c
      h := h * delta
      if (delta - 1).abs < 2e-15 then return clamp01 (← finite (1 - factor * h))
  throw "Incomplete-gamma iteration budget exhausted."

private def betaFraction (a b x : Float) : Result Float := do
  let qab := a + b
  let qap := a + 1
  let qam := a - 1
  let mut c := 1.0
  let mut d := 1 - qab * x / qap
  if d.abs < 1e-300 then d := 1e-300
  d := 1 / d
  let mut h := d
  for m in [1:10001] do
    let fm := nf m
    let m2 := 2 * fm
    let mut aa := fm * (b - fm) * x / ((qam + m2) * (a + m2))
    d := 1 + aa * d
    if d.abs < 1e-300 then d := 1e-300
    c := 1 + aa / c
    if c.abs < 1e-300 then c := 1e-300
    d := 1 / d
    h := h * d * c
    aa := -(a + fm) * (qab + fm) * x / ((a + m2) * (qap + m2))
    d := 1 + aa * d
    if d.abs < 1e-300 then d := 1e-300
    c := 1 + aa / c
    if c.abs < 1e-300 then c := 1e-300
    d := 1 / d
    let delta := d * c
    h := h * delta
    if (delta - 1).abs < 2e-15 then return (← finite h)
  throw "Incomplete-beta iteration budget exhausted."

def betaCDF (a b x : Float) : Result Float := do
  if !a.isFinite || !b.isFinite || a <= 0 || b <= 0 || !x.isFinite then
    throw "Invalid incomplete-beta arguments."
  if x <= 0 then return 0
  if x >= 1 then return 1
  if a > 100000 || b > 100000 then throw "Incomplete-beta shapes exceed the numerical limit (100000)."
  let front := ((← logGamma (a+b)) - (← logGamma a) - (← logGamma b)
    + a * x.log + b * (1-x).log).exp
  if x < (a+1)/(a+b+2) then
    return clamp01 (← finite (front * (← betaFraction a b x) / a))
  return clamp01 (← finite (1 - front * (← betaFraction b a (1-x)) / b))

/-- Cody's normal-tail approximation (same coefficients as core/src/stats/normal.cpp). -/
def normalCDF (x : Float) : Float := Id.run do
  if x.isNaN then return x
  let y := x.abs
  if y <= 0.67448975 then
    let a : Array Float := #[2.2352520354606839287, 161.02823106855587881,
      1067.6894854603709582, 18154.981253343561249, 0.065682337918207449113]
    let b : Array Float := #[47.20258190468824187, 976.09855173777669322,
      10260.932208618978205, 45507.789335026729956]
    let z := x*x
    let mut num := a[4]! * z
    let mut den := z
    for i in [:3] do
      num := (num + a[i]!) * z
      den := (den + b[i]!) * z
    return 0.5 + x * (num + a[3]!) / (den + b[3]!)
  let mut tail := 0.0
  if y <= 5.6568542494923801952 then
    let c : Array Float := #[0.39894151208813466764, 8.8831497943883759412,
      93.506656132177855979, 597.27027639480026226, 2494.5375852903726711,
      6848.1904505362823326, 11602.651437647350124, 9842.7148383839780218,
      1.0765576773720192317e-8]
    let d : Array Float := #[22.266688044328115691, 235.38790178262499861,
      1519.377599407554805, 6485.558298266760755, 18615.571640885098091,
      34900.952721145977266, 38912.003286093271411, 19685.429676859990727]
    let mut num := c[8]! * y
    let mut den := y
    for i in [:7] do
      num := (num + c[i]!) * y
      den := (den + d[i]!) * y
    tail := (-0.5*y*y).exp * (num + c[7]!) / (den + d[7]!)
  else if y < 39 then
    let p : Array Float := #[0.21589853405795699, 0.1274011611602473639,
      0.022235277870649807, 0.001421619193227893466, 2.9112874951168792e-5,
      0.02307344176494017303]
    let q : Array Float := #[1.28426009614491121, 0.468238212480865118,
      0.0659881378689285515, 0.00378239633202758244, 7.29751555083966205e-5]
    let z := 1/(y*y)
    let mut num := p[5]! * z
    let mut den := z
    for i in [:4] do
      num := (num + p[i]!) * z
      den := (den + q[i]!) * z
    tail := (-0.5*y*y).exp * (0.39894228040143267794 - z*(num+p[4]!)/(den+q[4]!)) / y
  return if x < 0 then tail else 1-tail

def tCDF (df x : Float) : Result Float := do
  if !df.isFinite || df <= 0 || !x.isFinite then throw "Invalid t-distribution arguments."
  let tail := 0.5 * (← betaCDF (df/2) 0.5 (df/(df+x*x)))
  return if x < 0 then tail else 1-tail

/-- Monotone inverse CDF on an explicitly supplied bracket. -/
def quantile (cdf : Float → Result Float) (p lo hi : Float) : Result Float := do
  if !p.isFinite || p <= 0 || p >= 1 || !lo.isFinite || !hi.isFinite || lo >= hi then
    throw "Quantile needs 0 < p < 1 and a finite increasing bracket."
  if (← cdf lo) > p || (← cdf hi) < p then throw "Quantile is outside the supplied bracket."
  let mut a := lo
  let mut b := hi
  for _ in [:100] do
    let m := a + (b-a)/2
    if (← cdf m) < p then a := m else b := m
  finite (a+(b-a)/2)

def normalQuantile (p : Float) : Result Float := do
  if p > 0.5 then return -(← quantile (fun x => finite (normalCDF x)) (1-p) (-40) 0)
  quantile (fun x => finite (normalCDF x)) p (-40) 0
def tQuantile (df p : Float) : Result Float := do
  if p > 0.5 then return -(← quantile (tCDF df) (1-p) (-1e8) 0)
  quantile (tCDF df) p (-1e8) 0

private def simpson (f : Float → Float) (a b : Float) (n : Nat) : Float := Id.run do
  let h := (b-a)/nf n
  let mut total := f a + f b
  for i in [1:n] do total := total + (if i % 2 == 0 then 2 else 4) * f (a + nf i * h)
  return total*h/3

/-- Noncentral t as a normal/chi-square mixture. Integration error is estimated
by refinement, not a certified bound. Inadequate convergence returns an error. -/
def noncentralTCDF (df delta x : Float) : Result Float := do
  if !df.isFinite || df <= 0 || !delta.isFinite || !x.isFinite then
    throw "Invalid noncentral-t arguments."
  if delta == 0 then return (← tCDF df x)
  let a := df/2
  let lg := (← logGamma a)
  let mut width := 8 / df.sqrt
  let mut bracketed := false
  for _ in [:32] do
    let below ← gammaCDF a (a * (-width).exp)
    let above ← gammaCDF a (a * width.exp)
    if below < 1e-11 && 1-above < 1e-11 then
      bracketed := true
      break
    width := width * 1.5
    if width > 680 then break
  if !bracketed then throw "Noncentral-t integration could not bound the omitted tails."
  let f := fun z =>
    let v := z.exp
    normalCDF (x * (z/2).exp - delta) * (a * (a.log + z) - a*v - lg).exp
  let mut previous := simpson f (-width) width 64
  let mut n := 128
  for _ in [:8] do
    let current := simpson f (-width) width n
    if n >= 256 && (current-previous).abs < 2e-9 then
      return clamp01 (← finite current)
    previous := current
    n := n*2
  throw "Noncentral-t quadrature did not converge."

end SSC.Numeric
