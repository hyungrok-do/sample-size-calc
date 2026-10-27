// Tests call the actual compiled C++ engine. There is no JavaScript simulator fallback.
import createEngine from '../engine/ssc_web.mjs';
import { parsePowerSettings } from '../power-settings.js';
const M = await createEngine();
export function native(request) {
  const text = JSON.stringify(request), size = M.lengthBytesUTF8(text) + 1, ptr = M._malloc(size);
  if (!ptr) throw new Error('WASM allocation failed');
  M.stringToUTF8(text, ptr, size);
  let out;
  try {
    out = M._ssc_dag_json(ptr, 0);
    if (!out) throw new Error('WASM did not return an answer');
    const value = JSON.parse(M.UTF8ToString(out));
    if (!value.ok) throw new Error(value.error);
    return value.result;
  } finally { if (out) M._ssc_free(out); M._free(ptr); }
}
export const simulateGraph = (graph, options = {}) => native({ op: 'simulate', graph, options });
export const analyzeReplicate = (graph, estimand, data, side = 'two-sided') => native({ op: 'analyze', graph, estimand, data, side });
export const populationEffect = (graph, estimand, options = {}) => native({ op: 'population', graph, estimand, options });
export const calibrateGraph = (graph, estimand, target, settings, seed) => native({ op: 'calibrate', graph, estimand, target, settings, seed });
export const noEffectGraph = (graph, estimand) => native({ op: 'nullGraph', graph, estimand });
export const runPower = (graph, estimand, settings) => native({ op: 'power', graph, estimand, settings: parsePowerSettings(settings) });
