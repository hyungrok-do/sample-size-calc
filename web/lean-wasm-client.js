import { extractCalculations } from './lean-calculation-results.js?v=dev';

export class BrowserLeanCompiler {
  constructor(progress) { this.progress = progress; this.worker = null; this.pending = new Map(); this.serial = 0; }
  async initialize() {
    if (this.initializing) return this.initializing;
    if (!crossOriginIsolated || typeof SharedArrayBuffer !== 'function') {
      if ('serviceWorker' in navigator) {
        await navigator.serviceWorker.register('sw.js');
        await navigator.serviceWorker.ready;
        if (!sessionStorage.getItem('ssc-lean-isolation-reload')) {
          sessionStorage.setItem('ssc-lean-isolation-reload', '1');
          this.progress('Preparing browser isolation. Reloading once…');
          location.reload();
          return new Promise(() => {});
        }
      }
      throw Error('This browser needs cross-origin isolation for Lean WebAssembly. Reload once, or use a browser that supports SharedArrayBuffer.');
    }
    sessionStorage.removeItem('ssc-lean-isolation-reload');
    this.initializing = new Promise((resolve, reject) => {
      const worker = this.worker = new Worker(new URL('./lean-wasm-worker.js?v=dev', import.meta.url));
      this.rejectInitialization = reject;
      this.initializationTimer = setTimeout(() => this.stop('Loading the browser compiler timed out. Retry to continue from cached files.'), 300000);
      worker.onerror = event => { if (this.worker === worker) this.stop(event.message || 'The browser could not start Lean WebAssembly.'); };
      worker.onmessage = ({ data }) => {
        if (this.worker !== worker) return;
        if (data.type === 'progress') this.progress(data.text);
        else if (data.type === 'ready') { clearTimeout(this.initializationTimer); this.rejectInitialization = null; this.ready = true; resolve(data); }
        else if (data.type === 'result' || data.type === 'error') {
          if (data.id == null || data.fatal) { this.stop(data.message); return; }
          const pending = this.pending.get(data.id);
          if (!pending) return;
          this.pending.delete(data.id); clearTimeout(pending.timer);
          if (data.type === 'error') pending.reject(Error(data.message)); else pending.resolve(extractCalculations(data.value));
        }
      };
      worker.postMessage({ type: 'initialize' });
    });
    return this.initializing;
  }
  async compile(source) {
    await this.initialize();
    return new Promise((resolve, reject) => {
      const id = ++this.serial;
      const timer = setTimeout(() => this.stop('Compilation exceeded its browser time limit.'), 120000);
      this.pending.set(id, { resolve, reject, timer });
      this.worker.postMessage({ type: 'compile', id, source });
    });
  }
  stop(message = 'Compilation stopped.') {
    clearTimeout(this.initializationTimer);
    this.rejectInitialization?.(Error(message)); this.rejectInitialization = null;
    this.worker?.terminate(); this.worker = null; this.ready = false; this.initializing = null;
    for (const pending of this.pending.values()) { clearTimeout(pending.timer); pending.reject(Error(message)); }
    this.pending.clear();
  }
}
