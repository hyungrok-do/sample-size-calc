// Short-lived startup helper. The parent authenticates both decoder assets and
// every compressed input before sending them here. Terminate this worker after
// installing the libraries so the decoder's retained WASM heap is released.
let decoder;
let requests = Promise.resolve();

async function handle(data) {
  if (data.type === 'initialize') {
    if (decoder) throw Error('The decompression worker is already initialized.');
    if (typeof data.moduleURL !== 'string' || !(data.wasmBinary instanceof Uint8Array))
      throw Error('Supply the authenticated decoder module and WASM bytes.');
    const moduleURL = new URL(data.moduleURL, self.location.href);
    if (moduleURL.origin !== self.location.origin)
      throw Error('The decompression module must have the same origin.');
    const module = await import(moduleURL.href);
    if (typeof module.default !== 'function' || typeof module.decompress !== 'function')
      throw Error('Unsupported decompression module.');
    // Pinned brotli-dec-wasm 2.3.2 takes explicit bytes through this object,
    // avoiding its default URL fetch and the deprecated positional API.
    await module.default({ module_or_path: data.wasmBinary });
    decoder = module.decompress;
    self.postMessage({ type: 'ready', id: data.id });
    return;
  }
  if (data.type !== 'decompress') throw Error('Unknown decompression request.');
  if (!decoder) throw Error('The decompression worker is not initialized.');
  if (!(data.bytes instanceof Uint8Array) || !Number.isSafeInteger(data.expectedBytes)
      || data.expectedBytes < 0 || data.expectedBytes > 512 * 1024 * 1024)
    throw Error('Invalid compressed asset or expected size.');
  const bytes = decoder(data.bytes);
  if (!(bytes instanceof Uint8Array) || bytes.byteLength !== data.expectedBytes)
    throw Error('Incorrect decompressed asset size.');
  self.postMessage({ type: 'result', id: data.id, bytes }, [bytes.buffer]);
}

self.onmessage = ({ data }) => {
  // Initialization imports asynchronously. Serialize requests so a queued
  // decode cannot observe a partially initialized module.
  requests = requests.then(() => handle(data)).catch(error => {
    self.postMessage({ type: 'error', id: data?.id,
      message: error?.message || 'Asset decompression failed.' });
  });
};
