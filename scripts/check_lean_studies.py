#!/usr/bin/env python3
"""Native Lean -> statistical IR -> compiled C++, against independent references."""
from fractions import Fraction
import json
import math
from pathlib import Path
import sys
import tempfile

import numpy as np
from scipy.stats import norm

from lean_equations import ROOT, compile_source
from check_lean_equations import cpp_run


def value(q):
    return float(Fraction(q))


def power(n, delta, variance, alpha):
    z = norm.isf(alpha/2)
    lam = abs(delta)/math.sqrt(variance)*math.sqrt(n)
    return norm.cdf(-z-lam)+norm.sf(z-lam)


def checked(source):
    r = compile_source(source)
    assert r['ok'], r.get('error')
    assert r['leanKernelChecked']
    for plan in r['results']:
        assert plan['kind'] == 'statistical_design'
        assert set(plan['axioms']) <= {'propext','Quot.sound','Classical.choice'}
    return r['results']


def verify_size(plan):
    stats, crit, out = plan['statistics'], plan['criterion'], plan['evaluation']
    delta, variance, alpha, target = value(stats['effect']), value(stats['variance']), value(crit['alpha']), value(crit['targetPower'])
    expected = next(n for n in range(crit['minN'], crit['maxN']+1) if power(n,delta,variance,alpha)>=target)
    assert out['selectedN'] == expected
    assert math.isclose(out['power'], power(expected,delta,variance,alpha), abs_tol=2e-12)
    assert math.isclose(out['criticalValue'], norm.isf(alpha/2), abs_tol=2e-10)
    if expected>crit['minN']:
        assert out['predecessorPower']<target
        assert math.isclose(out['predecessorPower'],power(expected-1,delta,variance,alpha),abs_tol=2e-12)
    assert out['kernelCertified'] is False
    return expected


def main():
    examples={}
    for path in sorted((ROOT/'examples/lean-studies').glob('*.lean')):
        plan=checked(path.read_text())[0]
        n=verify_size(plan)
        examples[path.stem]=plan
        standalone=cpp_run(plan)
        assert standalone == plan['evaluation']
        print(path.stem, n, plan['sizeUnit'], 'power=',plan['evaluation']['power'], flush=True)
    repeated=examples['repeated-measures']
    assert value(repeated['statistics']['variance']) == 1.6
    assert repeated['evaluation']['selectedN']==51 and repeated['evaluation']['totalParticipants']==102
    assert examples['cluster-trial']['evaluation']['selectedN']==18
    assert examples['cluster-trial']['evaluation']['totalParticipants']==720
    bad_icc=(ROOT/'examples/lean-studies/cluster-trial.lean').read_text().replace('def icc : ℚ := 1/20', 'def icc : ℚ := 2')
    assert not compile_source(bad_icc)['ok']
    assert value(examples['targeted-ate']['statistics']['effect'])==0.12
    assert value(examples['targeted-ate']['statistics']['variance'])==0.6906
    # Independent Monte Carlo check of the multivariate outcome model, not of
    # the exporter. Vary correlation: ignoring it should change the sample size.
    rng=np.random.default_rng(20261103)
    model=repeated['model']; n=repeated['evaluation']['selectedN']
    covariance=np.array([[value(x) for x in row] for row in model['covariance']])
    c=np.array([value(x) for x in model['contrast']]); mean=np.array([value(x) for x in model['means']])
    means=rng.multivariate_normal(mean,covariance/n,size=50000)
    z=(means@c)/math.sqrt((c@covariance@c)/n)
    empirical=np.mean(np.abs(z)>norm.isf(.025))
    theoretical=repeated['evaluation']['power']
    assert abs(empirical-theoretical) < 5*math.sqrt(theoretical*(1-theoretical)/len(z))
    source=(ROOT/'examples/lean-studies/repeated-measures.lean').read_text()
    independent=checked(source.replace('3/5','0'))[0]
    assert verify_size(independent)>n
    # Multiple actual Lean commands, reuse checked arithmetic in a later proof.
    module=source+'\nsolve power given_size using study at 100\nsolve effect detectable using study at 100\nexample : variance study = 8/5 := repeated_measures_variance\n'
    plans=checked(module)
    assert len(plans)==3
    assert math.isclose(plans[1]['evaluation']['power'],power(100,.5,1.6,.05),abs_tol=2e-12)
    effect=plans[2]['evaluation']['detectableMagnitude']
    assert math.isclose(power(100,effect,1.6,.05),.8,abs_tol=2e-12)
    # An unseen three-dimensional contrast, generated from an SPD covariance.
    arbitrary='''open DesignCompiler.Planning
    def power (n : ℚ) : ℚ := n+1
    def effect : ℚ := 3/4
    def makeCov (a : ℚ) : List (List ℚ) := [[2,a,0],[a,3,1/2],[0,1/2,1]]
    def study : Study := {
      name := "Unseen contrast", estimand := "mu1-2*mu2+mu3", analysis := .knownCovarianceZ
      model := .gaussianContrast {
        labels := ["a","b","c"]
        means := [0,0,effect]
        coefficients := [1,-2,1]
        covariance := makeCov (1/2)
      }
      criterion := { alpha := 1/100, targetPower := 9/10, minN := 3, maxN := 10000 }
    }
    solve sample_size custom using study
    '''
    custom=checked(arbitrary)[0]
    assert value(custom['statistics']['variance'])==11
    verify_size(custom)
    # Singular valid covariance must remain usable for a positive-variance contrast.
    singular=arbitrary.replace('[[2,a,0],[a,3,1/2],[0,1/2,1]]','[[1,1,1],[1,1,1],[1,1,1]]').replace('[1,-2,1]', '[1,1,1]')
    assert value(checked(singular)[0]['statistics']['variance'])==9
    # Preserve the standard-error scale when V/n would underflow binary64.
    tiny = checked('''open DesignCompiler.Planning
def study : Study := {
  name := "Small effect scale", estimand := "Scaled estimand", analysis := .wald
  model := .asymptoticWald {
    alternative := 1/(10^162)
    variance := 1/(10^323)
    varianceSource := "Rescaled influence distribution"
    assumptions := ["Independent sampling and asymptotic linearity"]
  }
}
solve sample_size small_scale using study
''')[0]
    tiny_n = verify_size(tiny)
    expected_se = math.sqrt(value(tiny['statistics']['variance']))/math.sqrt(tiny_n)
    assert math.isclose(tiny['evaluation']['standardError'], expected_se, rel_tol=1e-12)
    # Invalid statistical specifications must not receive a power result.
    invalid=[
      (source.replace('[1, 3/5, 0, 0]', '[1, 2, 0, 0]').replace('[3/5, 1, 0, 0]','[2, 1, 0, 0]'), 'positive semidefinite'),
      (source.replace('[1, 3/5, 0, 0]', '[1, 1/2, 0, 0]'), 'symmetric'),
      (source.replace('coefficients := [1, -1, -1, 1]', 'coefficients := [1, -1]'), 'must agree'),
      (source.replace('analysis := .knownCovarianceZ', 'analysis := .tmle'), 'requires knownCovarianceZ'),
      (source.replace('targetPower := 4/5', 'targetPower := 1/100'), 'alpha'),
      (source.replace('minN := 2', 'minN := 0'), 'integer search'),
      ((ROOT/'examples/lean-studies/targeted-ate.lean').read_text().replace('propensity := 1/2','propensity := 0'), 'positivity'),
      ((ROOT/'examples/lean-studies/targeted-ate.lean').read_text().replace('weight := 3/5','weight := 1/2'), 'sum exactly'),
      ((ROOT/'examples/lean-studies/targeted-ate.lean').read_text().replace('solve sample_size targeted_ate using study','solve effect targeted_ate using study at 100'), 'change the planning strata'),
    ]
    for text,message in invalid:
        r=compile_source(text)
        assert not r['ok'] and message in r['error'], r
    no_effect=checked(source.replace('means := [0, 0, 0, 1/2]', 'means := [0, 0, 0, 0]'))[0]
    assert no_effect['evaluation']['status']=='not_found_in_domain'
    too_small=checked(source.replace('maxN := 10000','maxN := 10'))[0]
    assert too_small['evaluation']['status']=='not_found_in_domain'
    # Timer applies to the entire native pipeline.
    assert not compile_source(source,timeout=0.01)['ok']
    print(f'Independent multivariate Gaussian simulation: empirical={empirical:.5f}, analytic={theoretical:.5f}')
    print('Statistical Lean: covariance, allocation, clusters, EIF decomposition, general contrasts, power/effect inversion, integer minimality, invalid-model rejection, exported C++, and simulation reference passed.')


if __name__=='__main__':
    main()
