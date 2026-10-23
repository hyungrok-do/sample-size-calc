// The page's handle on the engine worker.  A cross-origin isolated page gets the engine that uses every core
// (ssc_web_mt); any other page, or a browser that cannot start it, gets the one-thread engine.  Stop ends the
// worker in the middle of a calculation and starts a fresh one.
const STARTUP_MS = 20000;      // from compiled to ready: starting the threads (a tenth of a second); more means never
const RECYCLE_BYTES = 64e6;    // an engine whose heap is this big after an answer is replaced once idle (it never shrinks)
const CRASHED = 'The engine stopped unexpectedly during this calculation, often for lack of memory with a very large '
  + 'answer. It has been restarted: try again, with fewer scenarios if this was a large one.';

export class Engine {
  constructor() {
    const asked = new URLSearchParams(location.search).get('engine');  // ?engine=single or ?engine=threads
    this.threaded = self.crossOriginIsolated === true && asked !== 'single';
    this.next = 1;
    this.failed = false;
    this.start();
  }

  start() {
    this.variant = this.threaded ? 'ssc_web_mt' : 'ssc_web';
    this.pending = new Map();
    this.recycle = false;
    this.failed = false;
    this.worker = new Worker(new URL(`./worker.js?engine=${this.variant}&v=dev`, import.meta.url), { type: 'module' });
    let watchdog = 0;
    this.worker.onmessage = ({ data }) => {
      if (data.stage === 'compiled') {
        if (this.threaded) watchdog = setTimeout(() => this.fallBack(), STARTUP_MS);
        return;
      }
      const call = this.pending.get(data.id);
      if (!call) return;
      this.pending.delete(data.id);
      if (data.crashed) {
        // A trapped engine cannot calculate again: a fresh one at once, before anyone can call the old one.  After
        // the multi-core engine, the one-thread engine for the rest of the visit, and the call may be run again.
        const retry = this.threaded;
        this.threaded = false;
        this.restart();
        call.reject(Object.assign(new Error(CRASHED), { crashed: true, retry, detail: data.error }));
        return;
      }
      if (data.ok) call.resolve(data);
      else call.reject(new Error(data.error));
      const v = data.value;
      if (v?.heapBytes !== undefined ? v.heapBytes > RECYCLE_BYTES : v?.memoryBytes > Math.max(RECYCLE_BYTES, v?.initialBytes ?? 0)) {
        this.recycle = true;
      }
      if (this.recycle && this.pending.size === 0) queueMicrotask(() => this.replace());
    };
    this.worker.onerror = (e) => {
      e.preventDefault();
      if (this.threaded) this.fallBack();
      else this.failAll(new Error(e.message || 'the engine could not start'));
    };
    this.ready = this.call('info').then(
      (r) => {
        clearTimeout(watchdog);
        return r.value;
      },
      (err) => {
        clearTimeout(watchdog);
        if (err.stopped || err.movedOn) return this.ready;  // a fresh worker took over: its answer
        if (this.threaded) {
          this.fallBack();
          return this.ready;
        }
        this.failed = true;
        throw err;
      },
    );
    this.onstart?.();
  }

  // The one-thread engine in place of a multi-core one that failed or never started.
  fallBack() {
    if (!this.threaded) return;
    this.threaded = false;
    this.restart();
  }

  // A fresh engine (offline too: the service worker keeps copies); the calls waiting for the old one go to it.
  restart() {
    const waiting = [...this.pending.values()];
    this.worker.terminate();
    const moved = new Error('moved to another engine');
    moved.movedOn = true;
    this.start();
    for (const call of waiting) {
      if (call.op === 'info') call.reject(moved);
      else this.post(call);
    }
  }

  // A fresh engine in place of an idle one that a large answer left big.  Not while offline: the new worker
  // could not download the engine.
  replace() {
    if (!this.recycle || this.pending.size || navigator.onLine === false) return;
    this.worker.terminate();
    this.start();
  }

  // After a failed start (such as a restart while offline): try again.
  retry() {
    if (!this.failed) return;
    this.worker.terminate();
    this.start();
  }

  failAll(err) {
    for (const call of this.pending.values()) call.reject(err);
    this.pending.clear();
  }

  post(call) {
    this.pending.set(call.id, call);
    this.worker.postMessage({ id: call.id, op: call.op, arg: call.arg });
  }

  call(op, arg) {
    return new Promise((resolve, reject) => this.post({ id: this.next++, op, arg, resolve, reject }));
  }

  async run(request, maxRows) {
    const { value } = await this.call('run', { request, maxRows });
    return value;
  }

  stop() {
    this.worker.terminate();
    const stopped = new Error('Stopped');
    stopped.stopped = true;
    this.failAll(stopped);
    this.start();
  }
}
