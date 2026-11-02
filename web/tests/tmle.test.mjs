import test from 'node:test';
import assert from 'node:assert/strict';
import { scenarioDesign, SCENARIOS } from '../dag-scenarios.js';
import { exampleGraph, parseGraph } from '../dag-model.js';
import { defaultTarget } from '../dag-power.js';
import { compilePowerPlan } from '../power-model.js';
import { exportDesign, parseDesign } from '../estimand-model.js';
import { parseAnalysis } from '../dag-functions.js';
import { native } from './dag-engine.mjs';
const near = (a, b, tol = 1e-8) => assert.ok(Math.abs(a-b)<tol, `${a} versus ${b}`);
const mean = (a) => a.reduce((s, v) => s+v, 0)/a.length;
const analyze = (graph, estimand, data, analysis) => native({ op: 'analyze', graph, estimand, data, analysis });
const simulate = (graph, n, seed = 42) => native({ op: 'simulate', graph, options: { n, seed } });

// Independent reference: saturated binary-W propensities are cell proportions;
// Q starts at the arm mean. Newton's method solves the logistic fluctuation,
// or the least-squares score gives its identity-link counterpart in closed form.
function reference(A, Y, pi, q0, q1, bounds = null, clusters = null) {
  const n=A.length, qs=[q0.slice(),q1.slice()], scaled=bounds ? Y.map((y)=>(y-bounds[0])/(bounds[1]-bounds[0])) : Y;
  for(let a=0;a<2;a++) {
    const h=pi.map((p)=>1/(a?p:1-p)), rows=A.flatMap((v,i)=>v===a?[i]:[]);
    let q=qs[a];
    if(bounds) {
      const offset=q.map((v)=>{ const p=Math.max(1e-6,Math.min(1-1e-6,(v-bounds[0])/(bounds[1]-bounds[0]))); return Math.log(p/(1-p)); });
      let epsilon=0;
      for(let k=0;k<100;k++) {
        q=offset.map((x,i)=>1/(1+Math.exp(-x-epsilon*h[i])));
        const score=rows.reduce((s,i)=>s+h[i]*(scaled[i]-q[i]),0), info=rows.reduce((s,i)=>s+h[i]**2*q[i]*(1-q[i]),0);
        if(Math.abs(score)<1e-11) break;
        epsilon+=score/info;
      }
      qs[a]=q.map((v)=>bounds[0]+(bounds[1]-bounds[0])*v);
    } else {
      const epsilon=rows.reduce((s,i)=>s+h[i]*(Y[i]-q[i]),0)/rows.reduce((s,i)=>s+h[i]**2,0);
      qs[a]=q.map((v,i)=>v+epsilon*h[i]);
    }
  }
  const estimate=mean(qs[1])-mean(qs[0]);
  const ic=Y.map((y,i)=>qs[1][i]-qs[0][i]+A[i]*(y-qs[1][i])/pi[i]-(1-A[i])*(y-qs[0][i])/(1-pi[i])-estimate);
  let variance;
  if(clusters) {
    const ids=[...new Set(clusters)], sums=ids.map((u)=>ic.reduce((s,v,i)=>s+(clusters[i]===u?v/n:0),0));
    variance=sums.reduce((s,v)=>s+v*v,0)*ids.length/(ids.length-1);
  } else variance=ic.reduce((s,v)=>s+(v-mean(ic))**2,0)/(n*(n-1));
  return { estimate, se:Math.sqrt(variance) };
}

test('continuous and bounded TMLE agree with an independent stratified reference and solve the EIF score', () => {
  const graph=parseGraph(exampleGraph('nonlinear-confounding'));
  const X=[],A=[],Y=[];
  for(const w of [0,1]) for(const a of [0,1]) for(let j=0;j<(w===a?18:6);j++) {
    X.push(w); A.push(a); Y.push(25+12*w+6*a+(j%3-1)*3);
  }
  const data={ n:A.length,columns:{X,A,Y} }, e=defaultTarget(graph), pi=X.map((w)=>w?.75:.25);
  const q=[0,1].map((a)=>A.map(()=>mean(Y.filter((_,i)=>A[i]===a))));
  for(const bounds of [null,[0,100]]) {
    const analysis={method:'tmle',basis:'custom',terms:[],propensity:{basis:'custom',terms:[{left:'X'}]},...(bounds?{outcomeBounds:bounds}:{})};
    const fit=analyze(graph,e,data,analysis), ref=reference(A,Y,pi,...q,bounds);
    near(fit.estimate,ref.estimate); near(fit.se,ref.se); assert.ok(fit.diagnostics.targetingScore<1e-10);
    assert.deepEqual(native({op:'plan',graph,estimand:e,analysis}),compilePowerPlan(graph,e,analysis));
  }
});

test('known bounds and explicit propensity truncation validate, persist and report their actual use', () => {
  const {graph,estimand,power}=scenarioDesign('tmle-bounded'), data=simulate(graph,300);
  const saved=parseDesign(exportDesign(graph,estimand,power)); assert.deepEqual(saved.power.analysis.outcomeBounds,[0,100]);
  const fit=analyze(graph,estimand,data,{method:'tmle',outcomeBounds:[0,100]}); assert.ok(fit.diagnostics.activeMean>0 && fit.diagnostics.activeMean<100);
  data.columns.Y[0]=101;
  assert.throws(()=>analyze(graph,estimand,data,{method:'tmle',outcomeBounds:[0,100]}),/bounds/);
  for(const value of [-.01,.5,NaN]) assert.throws(()=>parseAnalysis({propensityBound:value}),/bound/);
  for(const b of [[1,0],[0],['0',1]]) assert.throws(()=>parseAnalysis({outcomeBounds:b}),/bounds/);
  const weak=scenarioDesign('weak-overlap');
  const result=analyze(weak.graph,weak.estimand,simulate(weak.graph,1000),weak.power.analysis);
  assert.ok(result.diagnostics.propensityTruncated>0); assert.ok(result.diagnostics.maxInverseWeight<=40+1e-9);
  assert.ok(result.diagnostics.targetingScore<1e-8);
});

test('clustered repeated TMLE equals a subject-weighted independent-unit reference including held-out cluster fits', () => {
  const {graph,estimand}=scenarioDesign('tmle-clustered'), d=simulate(graph,48,12), T=graph.study.times.length;
  const A=d.columns.A.filter((_,i)=>i%T===0), ids=d.columns.cluster_id.filter((_,i)=>i%T===0);
  const Y=A.map((_,i)=>mean(d.columns.Y.slice(i*T,(i+1)*T)));
  for(const folds of [1,3]) {
    const ranks=[0,0], uf=new Map();
    ids.forEach((u,i)=>{ if(!uf.has(u)) uf.set(u,ranks[A[i]]++%folds); });
    const pi=[],q0=[],q1=[];
    for(let i=0;i<A.length;i++) {
      const train=A.map((_,j)=>j).filter((j)=>folds===1||uf.get(ids[i])!==uf.get(ids[j]));
      pi.push(mean(train.map((j)=>A[j])));
      q0.push(mean(train.filter((j)=>A[j]===0).map((j)=>Y[j]))); q1.push(mean(train.filter((j)=>A[j]===1).map((j)=>Y[j])));
    }
    const analysis={method:'tmle',folds,basis:'custom',terms:[],propensity:{basis:'custom',terms:[]}};
    const fit=analyze(graph,estimand,d,analysis), ref=reference(A,Y,pi,q0,q1,null,ids);
    near(fit.estimate,ref.estimate); near(fit.se,ref.se); assert.equal(fit.df,46);
    assert.equal(fit.diagnostics.independentUnits,48); assert.equal(fit.diagnostics.subjects,A.length);
  }
});

test('binary repeated summaries use derived bounds for signed contrasts and never read missing outcomes', () => {
  const {graph,estimand}=scenarioDesign('tmle-repeated-binary'), d=simulate(graph,500);
  for(const repeated of [{contrast:'average'},{contrast:'change'},{contrast:'custom',weights:[-2,0,1,3]}]) {
    const e={...estimand,repeated}, plan=native({op:'plan',graph,estimand:e,analysis:{method:'tmle'}});
    assert.deepEqual(plan.outcomeBounds,repeated.contrast==='average'?[0,1]:repeated.contrast==='change'?[-1,1]:[-2,4]);
    const fit=analyze(graph,e,d,{method:'tmle',folds:3}); assert.ok(fit.diagnostics.targetingScore<1e-8);
    assert.ok(fit.diagnostics.activeMean>=plan.outcomeBounds[0] && fit.diagnostics.activeMean<=plan.outcomeBounds[1]);
  }
  d.columns.Y_observed[3]=0; d.columns.Y[3]=999999;
  assert.throws(()=>analyze(graph,estimand,d,{method:'tmle'}),/observed|recorded|missing/i);
  graph.study.observation.mode='dropout';
  for(const plan of [()=>compilePowerPlan(graph,estimand,{method:'tmle'}),()=>native({op:'plan',graph,estimand,analysis:{method:'tmle'}})]) assert.throws(plan,/complete outcomes/);
});

test('TMLE rejects post-treatment adjustment, time-varying treatment and unrelated estimand families', () => {
  const {graph,estimand}=scenarioDesign('tmle-repeated');
  const bad={method:'tmle',basis:'custom',terms:[{left:'Time'}]};
  assert.throws(()=>compilePowerPlan(graph,estimand,bad),/before treatment/);
  assert.throws(()=>native({op:'plan',graph,estimand,analysis:bad}),/before treatment/);
  const g=structuredClone(graph); g.nodes.find((n)=>n.id==='A').scope='visit';
  assert.throws(()=>native({op:'plan',graph:g,estimand,analysis:{method:'tmle'}}),/baseline/);
  const survival=scenarioDesign('survival'); assert.throws(()=>native({op:'plan',graph:survival.graph,estimand:survival.estimand,analysis:{method:'tmle'}}),/nonsurvival/);
});

test('continuous TMLE recovers with either nuisance model and exposes double misspecification', () => {
  const {graph,estimand,power}=scenarioDesign('tmle-continuous');
  const run=(basis,propensity)=>native({op:'power',graph,estimand,settings:{...power,sampleSizes:[800],repetitions:180,populationSize:3000,seed:90,analysis:{method:'tmle',folds:3,basis,propensity:{basis}},checkNull:true,...(propensity?{analysis:{method:'tmle',folds:3,basis,propensity:{basis:propensity}}}:{})}});
  const q=run('quadratic','linear'),g=run('linear','quadratic'),both=run('linear','linear');
  for(const r of [q,g]) { const s=r.results[0]; assert.ok(Math.abs(s.alternative.bias)<.08); assert.ok(s.alternative.coverage>.88); assert.ok(s.null.probability<.12); }
  assert.ok(Math.abs(both.results[0].alternative.bias)>.1);
});

test('every catalog entry carries a compatible, explicit design, and new scenarios execute in C++', () => {
  assert.equal(SCENARIOS.length,28);
  for(const s of SCENARIOS) {
    const {graph,estimand,power}=scenarioDesign(s.id);
    if(!estimand) { assert.equal(s.id,'empty'); continue; }
    const nativePlan=native({op:'plan',graph,estimand,analysis:power.analysis});
    assert.equal(nativePlan.method,compilePowerPlan(graph,estimand,power.analysis).method,s.id);
    const fit=analyze(graph,estimand,simulate(graph,graph.study?.clusterSizes.length?48:s.id==='rare-binary'?2000:500,17),power.analysis);
    assert.ok(Number.isFinite(fit.se) && fit.se>0,s.id);
  }
});

test('TMLE power diagnostics aggregate deterministically across threads and retain failed rare-event fits', () => {
  const {graph,estimand,power}=scenarioDesign('tmle-clustered');
  const req={op:'power',graph,estimand,settings:{...power,sampleSizes:[30],repetitions:30,populationSize:1000}};
  const one=native({...req,threads:1}),multi=native({...req,threads:3});
  assert.deepEqual(one.results,multi.results);
  assert.ok(one.results[0].alternative.targetedDiagnostics.maxInverseWeight>1);
  const rare=scenarioDesign('rare-binary'); rare.graph.nodes.find((n)=>n.id==='Y').intercept=-8;
  const r=native({op:'power',graph:rare.graph,estimand:rare.estimand,settings:{...rare.power,sampleSizes:[50],repetitions:20,populationSize:1000}}).results[0].alternative;
  assert.ok(r.failures>0); assert.equal(r.repetitions,20); assert.equal(r.valid+r.failures,20); assert.equal(r.probability,r.rejected/20);
});
